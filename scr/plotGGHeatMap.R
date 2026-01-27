plotGGHeatMap <- function(g.list, mm.counts, dset, pway.name){
  
  #' function to plot and a save a heatmap of scaled miroarray intensity using GGPlot2
  #'
  #' @param g.list list of genes (in the pathway) to filter the matrix on - first variable to make mapping easy
  #' @param mm.counts matrix of scale intensities
  #' @param dset the dataset name, for saving in the appropriate folder
  #' @param pway.name the name of the pathway being plotted for file naming 
  
  # filter on the gene list, and add sample display labels - dataframe is still grouped
  counts.mat.gene <- mm.counts %>%
    dplyr::select(sampleID, c(which(colnames(.) %in% g.list), condition)) %>% 
    dplyr::arrange(c("condition")) %>% 
    dplyr::mutate(prefix = case_when(condition == 0 ~ "c",
                                     condition == 1 ~ "s",
                                     T ~ NA),
                  id = row_number(),
                  display.id = paste(prefix, id, sep="")) %>%
    dplyr::select(-c(prefix,id))
  
  # make anon sample ID an ordered factor
  counts.mat.gene$display.id <- factor(counts.mat.gene$display.id, levels = counts.mat.gene$display.id)
  

  # pivot longer and plot heat map
  
  heatmap <- counts.mat.gene %>% 
    tidyr::pivot_longer(cols = -c(sampleID, condition, display.id), names_to = 'gene.name', values_to = 'exp') %>% 
    ggplot(aes(x = display.id,  y = gene.name, fill = exp)) + 
      geom_tile() +
      scale_fill_gradientn(limits=c(0, 1), breaks = seq(-1,1,by=0.20), colours = colorRampPalette(brewer.pal(9, "BuGn"))(20)) +
      coord_equal() +
      theme(axis.text.x = element_text(angle = 270)) +
      xlab("Patient") +
      ylab("Gene") +
      labs(fill = element_blank()) +
      theme(
            # legend.position = 'none',
            legend.position = 'right',
            axis.text.x = element_blank(),
            axis.text.y = element_text(size = 20),
            axis.title.x = element_text(size = 24),
            axis.title.y = element_text(size = 24), 
            legend.text = element_text(size = 16))
  
  # print(heatmap)
   ggsave(filename = here("figs", "pathway_het", paste(dset, pway.name, "ggheatmap", "pdf", sep = ".")),
         plot = heatmap,
         width = 12,
         height = 4)
  
}
