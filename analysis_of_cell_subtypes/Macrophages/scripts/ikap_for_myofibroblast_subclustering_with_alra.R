BottomUpMerge <- function (sobj, k.max, npc, random.seed) 
{
  k.clustering <- 0
  clusterings <- list()
  clust.r <- 1
  sobj <- FindClusters(object = sobj, reduction.type = "pca", 
                       dims.use = 1:npc, resolution = clust.r, print.output = 0, 
                       save.SNN = TRUE, random.seed = random.seed)
  cat("Iteration for nPC =", npc, ", r = 1.0")
  while (length(unique(sobj@active.ident)) < k.max) {
    clust.r <- clust.r + 0.2
    cat(",", clust.r)
    sobj <- FindClusters(object = sobj, reduction.type = "pca", 
                         dims.use = 1:npc, resolution = clust.r, print.output = 0, 
                         save.SNN = TRUE, random.seed = random.seed)
  }
  cat("\n")
  clusterings[[(k.clustering <- length(unique(sobj@active.ident)))]] <- as.character(sobj@active.ident)
  while (k.clustering > 2) {
    merged <- NearestCluster(sobj@reductions$pca@cell.embeddings[, 
                                                                 1:npc], clusterings[[k.clustering]])
    clustering.merged <- clusterings[[k.clustering]]
    clustering.merged[which(clustering.merged == merged[1])] <- merged[2]
    clusterings[[k.clustering - 1]] <- as.character(as.integer(as.factor(clustering.merged)))
    k.clustering <- k.clustering - 1
  }
  clusterings[[1]] <- rep("1", nrow(sobj@reductions$pca@cell.embeddings))
  return(clusterings[1:k.max])
}


ComputeMarkers <- function (sobj, gap.gain, candidates, out.dir) 
{
  markers.all <- list()
  out.xls <- list()
  out.xls$gap.gain <- cbind(PC_K = rownames(gap.gain), gap.gain)
  for (i in 1:length(candidates$k)) {
    clustering.label <- paste0("PC", candidates$pc[i], 
                               "K", candidates$k[i])
    Idents(sobj) <- clustering.label
    sobj <- RunTSNE(sobj, dims.use = 1:candidates$pc[i])
    ggsave(TSNEPlot(sobj), filename = paste0(out.dir, 
                                             "/", clustering.label, "_tSNE.pdf"))
    sobj.markers <- FindAllMarkers(object = sobj, only.pos = TRUE, 
                                   min.pct = 0.25, thresh.use = 0.25)
    sobj.markers$AUROC <- NA
    for (j in 1:nrow(sobj.markers)) {
      sobj.markers$AUROC[j] <- roc.curve(scores.class0 = sobj@assays$alra@data[sobj.markers$gene[j], 
      ], weights.class0 = sobj@active.ident == sobj.markers$cluster[j])$auc
    }
    top.10 <- sobj.markers %>% group_by(cluster) %>% top_n(10, 
                                                           avg_log2FC)
    ggsave(DoHeatmap(object = sobj, features = top.10$gene), 
           filename = paste0(out.dir, "/", clustering.label, 
                             "_DE_genes_LCF.pdf"), units = "in", 
           width = 12, height = 8)
    out.xls[[clustering.label]] <- sobj.markers[, c("gene", 
                                                    "p_val", "avg_log2FC", "pct.1", 
                                                    "pct.2", "p_val_adj", "cluster", 
                                                    "AUROC")]
    markers.all[[clustering.label]] <- sobj.markers
  }
  write.xlsx(out.xls, file  = paste0(out.dir, "/data.xlsx"))
  saveRDS(markers.all, file = paste0(out.dir, "/markers.all.rds"))
  return(markers.all)
}





DecisionTree <- function (sobj, markers, out.dir, plot.decision.tree) 
{
  data.rpart <- list()
  for (candidate in names(markers)) {
    data.rpart[[candidate]] <- list()
    genes.candidate <- unique(markers[[candidate]]$gene[which(markers[[candidate]]$p_val_adj < 
                                                                0.01)])
    for (clust in as.character(unique(markers[[candidate]]$cluster))) {
      if (length(genes.candidate) == 0) {
        next
      }
      else if (length(genes.candidate) == 1) {
        data <- data.frame(as.factor(sobj@meta.data[[candidate]] == 
                                       clust), as.numeric(sobj@assays$alra@data[genes.candidate, 
                                       ]))
        colnames(data) <- c("label", genes.candidate)
      }
      else {
        data <- as.data.frame(t(as.matrix(sobj@assays$alra@data[genes.candidate, 
        ])))
        data$label <- as.factor(sobj@meta.data[[candidate]] == 
                                  clust)
      }
      data.rpart[[candidate]][[clust]] <- rpart(label ~ 
                                                  ., data = data)
    }
  }
  summary.rpart <- data.frame(matrix(NA, nrow = length(markers), 
                                     ncol = 20))
  rownames(summary.rpart) <- names(markers)
  colnames(summary.rpart) <- paste0("s", 1:20)
  for (candidate in names(data.rpart)) {
    for (nsplit in 1:20) {
      err.rate <- 0
      for (clust in unique(as.character(sobj@meta.data[[candidate]]))) {
        if (is.null(data.rpart[[candidate]][[clust]])) {
          err.rate <- err.rate + 1
        }
        else {
          err.rate <- err.rate + data.rpart[[candidate]][[clust]]$cptable[max(which(data.rpart[[candidate]][[clust]]$cptable[, 
                                                                                                                             "nsplit"] <= nsplit)), "rel error"] * 
            length(which(sobj@meta.data[[candidate]] == 
                           clust))/nrow(sobj@meta.data)
        }
      }
      summary.rpart[candidate, nsplit] <- err.rate
    }
  }
  pdf(paste0(out.dir, "/DT_plot.pdf"))
  if (plot.decision.tree) {
    for (candidate in names(data.rpart)) {
      clust.sorted <- unique(as.character(sobj@meta.data[[candidate]]))
      if (!any(is.na(suppressWarnings(as.integer(clust.sorted))))) 
        clust.sorted <- as.character(sort(as.integer(clust.sorted)))
      for (clust in clust.sorted) {
        if (!is.null(data.rpart[[candidate]][[clust]])) {
          rpart.plot(data.rpart[[candidate]][[clust]], 
                     roundint = F, main = paste0(candidate, ":", 
                                                 clust))
        }
      }
    }
  }
  dev.off()
  saveRDS(data.rpart, file = paste0(out.dir, "/DT.rds"))
  saveRDS(summary.rpart, file = paste0(out.dir, "/DT_summary.rds"))
  return(summary.rpart)
}




ikap3 <- function (sobj, pcs = NA, pc.range = 20, k.max = NA, r.kmax.est = 1.5, 
                   out.dir = "./IKAP", scale.data = F, confounders = c("mitoRatio", "CC.Difference"),
                   plot.decision.tree = TRUE, random.seed = 0) 
{
  dir.create(out.dir, recursive = T)
  if (scale.data) {
    if (!all(confounders %in% colnames(sobj@meta.data))) {
      warning(confounders[which(!confounders %in% colnames(sobj@meta.data))], 
              "not in Seurat metadata: skipped for regression.\n")
    }
    confounders <- intersect(confounders, colnames(sobj@meta.data))
    if (length(confounders) > 0) 
      sobj <- ScaleData(sobj, vars.to.regress = confounders)
    #sobj <- SCTransform(sobj, vars.to.regress = confounders)
    else sobj <- ScaleData(sobj)
  }
  cat("Finding variable genes for clustering ... \n")
  # sobj <- FindVariableFeatures(sobj)
  cat("Running PCA ... \n")
  if (is.na(pcs)) {
    sobj <- RunPCA(sobj, pcs.compute = 50, do.print = F)
    pc.change <- which(abs(diff(sobj@reductions$pca@stdev)/sobj@reductions$pca@stdev[2:length(sobj@reductions$pca@stdev)]) > 
                         0.1)
    while (length(pc.change) > 0 && max(pc.change) + pc.range + 
           2 > length(sobj@reductions$pca@stdev)) {
      sobj <- RunPCA(sobj, pcs.compute = max(pc.change) + 
                       pc.range + 2 + 10, do.print = F)
      pc.change <- which(abs(diff(sobj@reductions$pca@stdev)/sobj@reductions$pca@stdev[2:length(sobj@reductions$pca@stdev)]) > 
                           0.1)
    }
    pcs <- if (length(pc.change) == 0) 
      2:(pc.range + 2)
    else (max(pc.change) + 2):(max(pc.change) + pc.range + 
                                 2)
  }
  else {
    sobj <- RunPCA(sobj, pcs.compute = max(pcs))
  }
  if (is.na(k.max)) {
    cat("Determine k.max.\n")
    sobj <- FindNeighbors(sobj, dims = 1:2, reduction="umap", force.recalc = T)
    sobj <- FindClusters(object = sobj, reduction.type = "pca", 
                         dims.use = 1:min(pcs), resolution = r.kmax.est, print.output = 0, 
                         save.SNN = TRUE, random.seed = random.seed)
    k.min.pc <- length(unique(sobj@active.ident))
    sobj <- FindClusters(object = sobj, reduction.type = "pca", 
                         dims.use = 1:max(pcs), resolution = r.kmax.est, print.output = 0, 
                         save.SNN = TRUE, random.seed = random.seed)
    k.max.pc <- length(unique(sobj@active.ident))
    k.max <- as.integer((k.min.pc + k.max.pc)/2)
    cat("k.max =", k.max, "\n")
  }
  gap.gain <- data.frame(matrix(NA, ncol = k.max - 1, nrow = length(pcs)))
  colnames(gap.gain) <- as.character(2:k.max)
  rownames(gap.gain) <- paste0(pcs)
  cat("Perform clustering for every nPC:\n")
  for (npc in pcs) {
    clusterings <- BottomUpMerge(sobj, k.max, npc, random.seed)
    gap.stat <- GapStatistic(sobj@reductions$pca@cell.embeddings[, 
                                                                 1:npc], clusterings)
    names(clusterings) <- paste0("PC", npc, "K", 
                                 1:k.max)
    sobj@meta.data <- cbind(sobj@meta.data, as.data.frame(clusterings)[, 
                                                                       2:k.max])
    gap.gain[as.character(npc), ] <- diff(gap.stat$gap)
  }
  candidates <- SelectCandidate(gap.gain)
  cat("Compute marker gene lists ... \n")
  markers.all <- ComputeMarkers(sobj, gap.gain, candidates, 
                                out.dir)
  cat("Build decision tree ... \n")
  summary.rpart <- DecisionTree(sobj, markers.all, out.dir, 
                                plot.decision.tree)
  cat("Plotting summary ... \n")
  PlotSummary(gap.gain, summary.rpart, markers.all, out.dir)
  return(sobj)
}
