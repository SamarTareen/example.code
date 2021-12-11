

# function to get a traditional volcano plot for DEGs
# copyright Samar H. K. Tareen

getVolcanoPlot <- function(FindMarkers_result, 
                           fold_change_thresh = 1, 
                           adj_p_val_thresh = 0.001){
  require(ggplot2)
  require(ggExtra)
  
  plot_data <- FindMarkers_result
  
  #finding which genes to colour
  plot_data$coloured <- "not_coloured"
  plot_data$coloured[which(
    plot_data$p_val_adj <= adj_p_val_thresh & 
      abs(plot_data$avg_log2FC) >= fold_change_thresh)] <- "coloured"
  plot_data <- plot_data[order(plot_data$coloured, decreasing = T),]
  
  #get legend labels and point colours
  if(length(unique(plot_data$coloured)) < 2){
    labels <- "Other genes"
    values <- "darkgrey"
  }
  else{
    labels <- c(paste0("Log |FC|\u2265 ", 
                       fold_change_thresh, 
                       "\nand adj.p.value\n\u2264 ", 
                       adj_p_val_thresh), 
                "Other genes")
    values <- c("#fd7239", "darkgrey")
  }
  
  #calculating x-axis so that inf values can be replaced
  x.axis <- plot_data$avg_log2FC
  x.axis.max <- x.axis
  x.axis.max[x.axis.max == Inf | x.axis.max == -Inf] <- NA
  x.axis.max <- max(abs(x.axis.max), na.rm = T)
  x.axis.max.nodetect <- x.axis.max + (0.5*x.axis.max)
  
  xbreaks <- 
    sort(c(-fold_change_thresh, fold_change_thresh,
           seq(0, (-x.axis.max)-1, -2),
           seq(0, (x.axis.max)+1, 2)))#, 
  #x.axis.max.nodetect, -x.axis.max.nodetect))
  if(x.axis.max.nodetect < max(max(xbreaks), abs(min(xbreaks)))){
    x.axis.max.nodetect <- max(max(xbreaks), abs(min(xbreaks))) + 
      (0.5*max(max(xbreaks), abs(min(xbreaks))))
  }
  xbreaks <- c(xbreaks, x.axis.max.nodetect, -x.axis.max.nodetect)
  xbreaks <- sort(unique(xbreaks))
  
  x.axis[x.axis == Inf] <- x.axis.max.nodetect
  x.axis[x.axis == -Inf] <- -x.axis.max.nodetect
  
  xlabels <- as.character(xbreaks)
  xlabels[xlabels == x.axis.max.nodetect | xlabels == -x.axis.max.nodetect] <- 
    sprintf("Detected\nonly in\nthis group")
  
  #calculating y-axis so that 0-p.values can be replaced by 
  # .Machine$double.xmin
  y.axis <- -log10(plot_data$p_val_adj)
  y.axis[y.axis == Inf] <- NA
  y.axis[is.na(y.axis)] <- max(y.axis, na.rm = T) + 
    (0.1*max(y.axis, na.rm = T))
  
  scatter_plot <- ggplot(data = plot_data,
                         aes(#x = plot_data$avg_logFC, 
                           x = x.axis, 
                           #y = -log10(plot_data$p_val_adj), 
                           y = y.axis, 
                           colour = plot_data$coloured, 
                           text = rownames(plot_data))) + 
    geom_hline(yintercept = 0, colour = "#4d4d4d", size = 0.25) + 
    geom_hline(yintercept = -log10(adj_p_val_thresh), 
               colour = "#4d4d4d", size = 0.25) + 
    geom_vline(xintercept = 0, colour = "#4d4d4d", size = 0.25) + 
    geom_vline(xintercept = fold_change_thresh, 
               colour = "#4d4d4d", size = 0.25) + 
    geom_vline(xintercept = -fold_change_thresh, 
               colour = "#4d4d4d", size = 0.25) + 
    geom_point(shape = 19, size = 1.5) + 
    scale_color_manual(
      labels = labels,
      values = values) +
    labs(colour = "Genes", 
         y = "-log10 of adjusted p.value", 
         x = "Average log2 fold change", 
         title = ""
    ) + 
    theme_light(base_size = 15) + 
    theme(plot.caption = element_text(hjust = 0, face= "italic"), 
          legend.key.width = unit(0.15,"line")
    )
  
  ybreaks <- sort(c(
    na.omit(ggplot_build(scatter_plot)$layout$panel_params[[1]]$y$breaks),
    -log10(adj_p_val_thresh)))
  # ybreaks <- sort(c(0, -log10(0.01)))
  
  ylabels <- as.character(ybreaks)
  ylabels[ylabels == -log10(adj_p_val_thresh)] <- paste0("\nadj.p\n", 
                                                         adj_p_val_thresh, 
                                                         "\n\n\n")
  #print(xbreaks)
  
  scatter_plot <- scatter_plot +
    scale_x_continuous(breaks = xbreaks, 
                       limits = c(min(xbreaks)-1, max(xbreaks)+1), 
                       labels = xlabels) +
    scale_y_continuous(breaks = ybreaks, labels = ylabels) + 
    expand_limits(x = c(-ceiling(max(abs(plot_data$avg_log2FC))), 
                        ceiling(max(abs(plot_data$avg_log2FC))))
    ) + 
    theme(panel.grid.minor.x = element_blank(), 
          #panel.grid.major.x = element_blank(), 
          panel.grid.minor.y = element_blank(), 
          #panel.grid.major.y = element_blank(),
          #axis.line = element_line(colour = 'black')
          panel.grid.major = element_line(colour = "#4d4d4d"),
          axis.ticks = element_line(colour = "#4d4d4d")
    )
  
  scatter_plot <- scatter_plot + ggExtra::removeGrid()
  
  return(scatter_plot)
  
}