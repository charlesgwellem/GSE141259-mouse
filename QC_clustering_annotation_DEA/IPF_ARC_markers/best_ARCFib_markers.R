library(Seurat)
library(ggplot2)

list.files()
seurat <- readRDS("subclustered_fibroblasts.rds")
DimPlot(seurat, group.by = "group")


DimPlot(seurat)


arc_markers <- FindMarkers(seurat,
                           only.pos = T, ident.1 = "stro_ARCFib")

arc_markers_sig <- arc_markers[arc_markers$p_val_adj < 0.05, ]


FeaturePlot(seurat, features = "SLC12A8", order = T,
            reduction = "umap")+
  theme(aspect.ratio = 1)

DotPlot(seurat, features = "SLC12A8")+
  theme(aspect.ratio = 1)

VlnPlot(seurat, features = "SLC12A8")+
  theme(aspect.ratio = 1)


ribos <- rownames(arc_markers_sig)[grep("RP[SL]",
                                        rownames(arc_markers_sig))]

arc_markers_sig <- setdiff(rownames(arc_markers_sig), ribos)

FeaturePlot(seurat, features = "TM4SF1", split.by = "group", order = T)

DefaultAssay(seurat) <- "SCT"

the_data <- FetchData(seurat, vars = c(arc_markers_sig, 
                                       "cellid_version3"))
head(the_data)

any(colnames(the_data)%in%"TM4SF1")

library(rpart)


gol.rp <- rpart(cellid_version3~., data=the_data, method="class", cp=0.001)
plot(gol.rp, branch=1,margin=0.1); text(gol.rp, digits=3, use.n=TRUE)


FeaturePlot(seurat, features = "MYO1E", order = T,
            reduction = "umap")+
  theme(aspect.ratio = 1)

vip::vip(gol.rp, 20)


df <- as.data.frame(gol.rp$variable.importance)
df[rownames(df)%in%"TM4SF1", ]

colnames(df)<-"varImp"

df$gene <- rownames(df)

df <- df[order(df$varImp, decreasing = T), ]

head(df, 20)


df_keep <- df[df$varImp>0.5, ]



DefaultAssay(seurat) <- "SCT"

the_data_floating <- FetchData(seurat, vars = c(df_keep$gene, 
                                       "cellid_version3"))

the_data_floating$cellid_version3 <- ifelse(
  the_data_floating$cellid_version3=="stro_ARCFib",
  1, 0
)

write.csv(the_data_floating, "arc_markers.csv",quote = F,
          row.names = F)

library(genefilter)

# find the most correlating genes to the most important separating gene
clostest_gene <- genefinder(as.matrix(GetAssayData(seurat, assay = "SCT")), 
                            c("ARC"), 
                            numResults=10, scale="none",
                            method="euclidean")


clostest_gene[[1]]$indices
round(clostest_gene[[1]]$dists,1)
near_genes <- rownames(seurat)[clostest_gene[[1]]$indices]
near_genes
FeaturePlot(seurat, features = c(near_genes)[1:4],
            order = T, reduction = "umap")

FeaturePlot(seurat, features = "HIPK2",
            order = T,  reduction = "umap")

d <- DotPlot(seurat, features = "UGDH")
d$data

RidgePlot(seurat, features = "TM4SF1")

best <- c('ARC', 'COL6A1', 'JARID2', 'IGF2', 'GPRC5A',
          'ABL2', 'LMNA', 'GLIS3', 'MAP4K4', 'ATP13A3',
          'DDX21', 'CRABP2', 'COL6A2', 'DKK1', 'TPT1',
          'ABLIM1', 'SNTG1', 'FKBP1C', 'IFI27', 'PCSK6')

VlnPlot(seurat, features = "SNTG1")



FeaturePlot(seurat, features = "SNTG1", reduction = "umap")+
  theme(aspect.ratio = 1)

h <- DotPlot(seurat, features = "SNTG1")
h$data
h + 
  theme (axis.text.x = element_text (angle = 45, vjust = 1, hjust=1)) 

# igraph
library(igraph)

DefaultAssay(seurat) <- "SCT"

# Filter the expression set object to include only genes of significant effect
highVarimp <- seurat[df_keep$gene, ]

matrix <- GetAssayData(highVarimp, assay = "SCT",
                       slot = "data")

matrix <- as.matrix(matrix)

library(igraph)

# Create a graph adjacency based on correlation distances between genes in  pairwise fashion.
g <- graph.adjacency(
  as.matrix(as.dist(cor(t(matrix), method="spearman"))),
  mode="undirected",
  weighted=TRUE,
  diag=FALSE
)


#g <- graph.adjacency(
#  as.matrix(dist(estrogenMainEffects, method="euclidean")),
#  mode="undirected",
#  weighted=TRUE,
#  diag=FALSE
#)

# Simplfy the adjacency object
g <- simplify(g, remove.multiple=TRUE, remove.loops=TRUE)

# Colour negative correlation edges as blue
E(g)[which(E(g)$weight<0)]$color <- "darkblue"

# Colour positive correlation edges as red
E(g)[which(E(g)$weight>0)]$color <- "darkred"

# Convert edge weights to absolute values
E(g)$weight <- abs(E(g)$weight)

# Change arrow size
# For directed graphs only
#E(g)$arrow.size <- 1.0

# Remove edges below absolute Pearson correlation 0.8
g <- delete_edges(g, E(g)[which(E(g)$weight<0.1)])

# Remove any vertices remaining that have no edges
g <- delete.vertices(g, degree(g)==0)

# Assign names to the graph vertices (optional)
V(g)$name <- V(g)$name

# Change shape of graph vertices
V(g)$shape <- "sphere"

# Change colour of graph vertices
V(g)$color <- "skyblue"

# Change colour of vertex frames
V(g)$vertex.frame.color <- "white"

# Scale the size of the vertices to be proportional to the level of expression of each gene represented by each vertex
# Multiply scaled vales by a factor of 10
scale01 <- function(x){(x-min(x))/(max(x)-min(x))}
vSizes <- (scale01(apply(matrix, 1, mean)) + 1.0) * 10

# Amplify or decrease the width of the edges
edgeweights <- E(g)$weight * 2.0

# Convert the graph adjacency object into a minimum spanning tree based on Prim's algorithm
mst <- mst(g, algorithm="prim")

# Plot the tree object
plot(
  mst,
  layout=layout.fruchterman.reingold,
  edge.curved=TRUE,
  vertex.size=vSizes,
  vertex.label.dist=-0.5,
  vertex.label.color="black",
  asp=FALSE,
  vertex.label.cex=0.6,
  edge.width=edgeweights,
  edge.arrow.mode=0,
  main="My first graph")



mst.communities <- edge.betweenness.community(mst, weights=NULL, directed=FALSE)
mst.clustering <- make_clusters(mst, membership=mst.communities$membership)
V(mst)$color <- mst.communities$membership + 1

par(mfrow=c(1,2))
plot(
  mst.clustering, mst,
  layout=layout.fruchterman.reingold,
  edge.curved=TRUE,
  vertex.size=vSizes,
  vertex.label.dist=-0.5,
  vertex.label.color="black",
  asp=FALSE,
  vertex.label.cex=0.6,
  edge.width=edgeweights,
  edge.arrow.mode=0,
  main="My first graph")

plot(
  mst,
  layout=layout.fruchterman.reingold,
  edge.curved=TRUE,
  vertex.size=vSizes,
  vertex.label.dist=-0.5,
  vertex.label.color="black",
  asp=FALSE,
  vertex.label.cex=0.6,
  edge.width=edgeweights,
  edge.arrow.mode=0,
  main="My first graph")

degs <- degree(g)

degs <- as.data.frame(degs)
degs$gene <- rownames(degs)
degs <- degs[order(degs$degs, decreasing = T), ]
head(degs)

write.xlsx(degs, "degrees_networks_arcfib_enriched_genes.xlsx",
           overwrite = T)

betw <- betweenness(g) 

betw <- as.data.frame(betw)
betw$gene <- rownames(betw)

betw <- betw[order(betw$betw, decreasing = T), ]
head(betw, 20)

write.xlsx(betw, "betweennes_networks_arcfib_enriched_genes.xlsx",
           overwrite = T)

library(openxlsx)


library(igraph)

# Calculate the out-degree of each vertex
g.outd <- degree(g, mode = c("out"))

# View a summary of out-degree
table(g.outd)

g.outd[order(g.outd, decreasing = T)]

# Make a histogram of out-degrees
hist(g.outd, breaks = 30)


library(igraph)

# Calculate betweenness of each vertex
g.b <- betweenness(g, directed = F)

# Show histogram of vertex betweenness
hist(g.b, breaks = 80)

# Create plot with vertex size determined by betweenness score
plot(g, 
     vertex.label = NA,
     edge.color = 'black',
     vertex.size = sqrt(g.b)+1,
     edge.arrow.size = 0.05,
     layout = layout_nicely(g))

# Find the vertex that has the maximum out-degree
which.max(g.outd)


# Make an ego graph
g184 <- make_ego_graph(g, diameter(g),
                       nodes = 'UGDH', mode = c("all"))[[1]]

# Get a vector of geodesic distances of all vertices from vertex 184 
dists <- distances(g184, "UGDH")

# Create a color palette of length equal to the maximal geodesic distance plus one.
colors <- c("black", "red", "orange", "blue", "dodgerblue", "cyan")

# Set color attribute to vertices of network g184.
V(g184)$color <- colors[dists+1]

# Visualize the network based on geodesic distance from vertex 184 (patient zero).
plot(g184, 
     vertex.label = dists, 
     vertex.label.color = "white",
     vertex.label.cex = .6,
     edge.color = 'black',
     vertex.size = 7,
     edge.arrow.size = .05,
     main = "Geodesic Distances from Patient Zero"
)




DotPlot(seurat, features = "HIPK2")+
  theme(aspect.ratio = 1)


RidgePlot(seurat, features = "HIPK2")+
  theme(aspect.ratio = 1)


VlnPlot(seurat, features = "HIPK2")+
  theme(aspect.ratio = 1)


FeaturePlot(seurat, features = "LRRC8C", split.by = "group",
            reduction = "umap", order = T)+
  theme(aspect.ratio = 1)


# Identify key nodes using eigenvector centrality
g.ec <- eigen_centrality(g)
which.max(g.ec$vector)

# Plot Forrest Gump Network
plot(g,
     vertex.label.color = "black", 
     vertex.label.cex = 0.6,
     vertex.size = 25*(g.ec$vector),
     edge.color = 'gray88',
     main = "Forrest Gump Network"
)


g.ec <- g.ec$vector[order(g.ec$vector, decreasing = T)]

g.ec <- as.data.frame(g.ec)
g.ec$gene <- rownames(g.ec)

write.xlsx(g.ec, "arc_enriched_genes_centralitiy.xlsx")

FeaturePlot(seurat, features = "HIPK2",
            reduction = "umap", order = T)+
  theme(aspect.ratio = 1)

RidgePlot(seurat, features = "DKK1")
