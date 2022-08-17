
select_spec_markers <- function(markers){
  
  specific_markers <- list()
  
  sig_markers <- subset(markers,
                        p_val_adj < 0.05 & avg_log2FC > 0)
  
  clusters <- unique(as.character(sig_markers$cluster))
  for (a_cluster in clusters) {
    cluster_oi <- subset(sig_markers,
                         cluster==a_cluster)
    other_clusters <- subset(sig_markers,
                             cluster!=a_cluster)
    cluster_oi_genes <- cluster_oi$gene
    other_clusters_genes <- other_clusters$gene
    specific_coi_genes <- setdiff(cluster_oi_genes,
                                  other_clusters_genes)
    specific_markers[[a_cluster]] <- specific_coi_genes
  }
  specific_markers
}
