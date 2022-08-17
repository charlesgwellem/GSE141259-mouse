cell_specific_ora <- function(marker_list,
                              dbs){
  
  for (cell_marker in names(marker_list)) {
    
    gene_list <- marker_list[[cell_marker]]
    
    ora_pways <- enrichr(gene_list, dbs_keep)
    
    
    cell_marker <- gsub("/", "", cell_marker)
    
    write.xlsx(ora_pways,
               paste0("results/", 
                      cell_marker, "_enriched_pways",
                      ".xlsx"),
               rowNames = TRUE)
    
    for (pway in names(ora_pways)) {
      
      pways_to_plot <- ora_pways[[pway]]
      
      pways_to_plot <- subset(pways_to_plot,
                              Adjusted.P.value < 0.05)
      
      if(nrow(pways_to_plot) != 0 & 
         !is.null(nrow(pways_to_plot))) {
        
        plotEnrich(pways_to_plot) +
          ggtitle(paste0(cell_marker, "_", pway)) +
          theme(aspect.ratio = 1) +
          scale_x_discrete(expand = c(0, 0)) +
          scale_y_continuous(expand = c(0, 0)) 
        ggsave(paste0("figures/", cell_marker, "_", pway,
                      ".pdf"),
               height = 6, width = 6, dpi = 500)
        
      } else {
        
        plotEnrich(pways_to_plot) +
          ggtitle(paste0(cell_marker, "_", pway)) +
          theme(aspect.ratio = 1)
        ggsave(paste0("figures/", cell_marker, "_", pway,
                      ".pdf"),
               height = 6, width = 6, dpi = 500)
        
      }
      
    }
    
  }
  
}