DimPlot(seurat, group.by = "manual_clusterwise_annot")

Idents(seurat) <- "manual_clusterwise_annot"

markers <- FindAllMarkers(seurat, only.pos = T)

dim(markers) # 37600     7

write.xlsx(markers, paste0(resultsdir, "/","manual_clusterwise_annot_res.res.2.4",
                           ".xlsx"), row.names = T)



sig_degs <- markers[markers$p_val_adj < 0.05, ]
dim(sig_degs) #34106     7

ribos <- sig_degs$gene[grep("Rp", sig_degs$gene)]

sig_genes <- setdiff(sig_degs$gene, ribos)

FeaturePlot(seurat, features = "Fn1", split.by = "group", order = T)

DefaultAssay(seurat) <- "SCT"

the_data <- FetchData(seurat, vars = c(sig_genes, 
                                       "manual_clusterwise_annot"))
head(the_data)

library(rpart)


gol.rp <- rpart(manual_clusterwise_annot~., data=the_data, method="class", cp=0.001)
plot(gol.rp, branch=1,margin=0.1); text(gol.rp, digits=3, use.n=TRUE)

FeaturePlot(seurat, features = "Ctsh", order = T)



# annotate again, using only the major cluster identities

DimPlot(seurat, group.by = "integrated_snn_res.0.2", label = T)

Idents(seurat) <- "integrated_snn_res.0.2"

markers_res.0.2 <- FindAllMarkers(seurat, only.pos = T)
dim(markers_res.0.2) # 15578     7



write.xlsx(markers_res.0.2, paste0(resultsdir, "/","markers_res.0.2",
                                   ".xlsx"), row.names = T)


sig_degs <- markers_res.0.2[markers_res.0.2$p_val_adj < 0.05, ]
dim(sig_degs) #14631     7

ribos <- sig_degs$gene[grep("Rp", sig_degs$gene)]

sig_genes <- setdiff(sig_degs$gene, ribos)
length(sig_genes) #6594

FeaturePlot(seurat, features = "Dntt", split.by = "group", order = T)

DefaultAssay(seurat) <- "SCT"

the_data <- FetchData(seurat, vars = c(sig_genes, 
                                       "integrated_snn_res.0.2"))
head(the_data)

dim(the_data) #18511  6595

library(rpart)


gol.rp <- rpart(integrated_snn_res.0.2~., data=the_data, method="class", cp=0.0001)
plot(gol.rp, branch=0,margin=0.1); text(gol.rp, digits=3, use.n=TRUE)

write.xlsx(gol.rp$frame, paste0(resultsdir, "/","gol.rp_integrated_snn_res.0.2",
                                ".xlsx"), row.names = T)


FeaturePlot(seurat, features = c("Ear2"), order = T)

ggsave(filename=paste0(figuresdir, "/", "Cotl1",
                       ".png"), width = 5, 
       height = 5, units = 'in', dpi = 300)


library(genefilter)

# find the most correlating genes to the most important separating gene
clostest_gene <- genefinder(as.matrix(GetAssayData(seurat, assay = "SCT")), 
                            c("Cox4i2"), 
                            numResults=25, scale="none", method="correlation")


clostest_gene[[1]]$indices
round(clostest_gene[[1]]$dists,1)
near_genes <- rownames(seurat)[clostest_gene[[1]]$indices]
near_genes
FeaturePlot(seurat, features = c(near_genes)[7:9],
            order = T)

FeaturePlot(bleo_saline, features = "Spp1",
            order = T, split.by = "group")

