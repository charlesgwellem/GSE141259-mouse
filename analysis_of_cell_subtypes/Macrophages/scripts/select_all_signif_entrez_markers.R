
select_all_sig_entrez_markers <- function(markers){
  
  all_sig_markers <- list()
  
  sig_markers <- subset(markers,
                        p_val_adj < 0.05 & avg_log2FC > 0)
  
  clusters <- unique(as.character(sig_markers$cluster))
  
  for (a_cluster in clusters) {
    
    cluster_oi <- subset(sig_markers,
                         cluster==a_cluster)
    
    other_clusters <- subset(sig_markers,
                             cluster!=a_cluster)
    
    cluster_oi_genes <- cluster_oi$ENTREZID
    
    other_clusters_genes <- other_clusters$ENTREZID
    
    all_sig_markers[[a_cluster]] <- cluster_oi_genes
    
  }
  
  all_sig_markers
}
