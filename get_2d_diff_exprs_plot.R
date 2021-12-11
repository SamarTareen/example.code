

# 2D differential expression plot (to combine GFP and IL2 volcano plots)
# copyright Samar H. K. Tareen

get_2d_diff_exprs_plot <- function(FindMarkers_result_x, 
                                   FindMarkers_result_y, 
                                   label_x = "", 
                                   label_y = "", 
                                   main_title = "", 
                                   label_quad1 = "", 
                                   label_quad3 = "", 
                                   #parse_labels = FALSE, 
                                   fold_change_thresh = 1, 
                                   adj_p_val_thresh = 0.001, 
                                   axisdh = 0.5, 
                                   show_corr = FALSE){
  require(ggplot2)
  require(ggExtra)
  
  #create the initial data structure for plotting
  plot_data <- merge(x = FindMarkers_result_x, y = FindMarkers_result_y, 
                     by = "row.names", all = T)
  rownames(plot_data) <- plot_data$Row.names
  
  #finding which genes to colour
  plot_data$coloured <- "not_coloured"
  plot_data$coloured[which(
    (plot_data$p_val_adj.x <= adj_p_val_thresh | 
       plot_data$p_val_adj.y <= adj_p_val_thresh) & 
      (abs(plot_data$avg_log2FC.x) >= fold_change_thresh | 
         abs(plot_data$avg_log2FC.y) >= fold_change_thresh))] <- "coloured"
  plot_data <- plot_data[order(plot_data$coloured, decreasing = T),]
  
  #get legend labels and point colours
  if(length(unique(plot_data$coloured)) < 2){
    labels <- "Other genes"
    values <- "darkgrey"
  }
  else{
    labels <- c(paste0("Log |FC|\u2265 ", 
                       fold_change_thresh, 
                       " and \nadj.p.value \u2264 ", 
                       adj_p_val_thresh, 
                       "\n(along either axis)"), 
                "Other genes")
    values <- c("#fd7239", "darkgrey")
  }
  
  #calculating x-axis so that inf values can be replaced
  x.axis <- c(plot_data$avg_log2FC.x, plot_data$avg_log2FC.y)
  x.axis.max <- x.axis
  x.axis.max[x.axis.max == Inf | x.axis.max == -Inf] <- NA
  x.axis.max <- max(abs(x.axis.max), na.rm = T)
  x.axis.max.nodetect <- x.axis.max + (0.2*x.axis.max)
  
  # # xbreaks <- 
  # #   sort(c(-fold_change_thresh, fold_change_thresh,
  # #          # seq(0, (-x.axis.max)-1, -2),
  # #          seq(0, (-x.axis.max)-1, -1),
  # #          # seq(0, (x.axis.max)+1, 2)))#, 
  # #          seq(0, (x.axis.max)+1, 1)))#, 
  # xbreaks <- 
  #   sort(c(-fold_change_thresh, fold_change_thresh,
  #          seq(0, floor(-x.axis.max/0.5)*0.5, -0.5),
  #          seq(0, ceiling(x.axis.max/0.5)*0.5, 0.5)))
  # xbreaks <-
  #   sort(c(-fold_change_thresh, fold_change_thresh,
  #          seq(0, floor(-x.axis.max/0.2)*0.2, -0.25),
  #          seq(0, ceiling(x.axis.max/0.2)*0.2, 0.25)))
  #           #x.axis.max.nodetect, -x.axis.max.nodetect))
  xbreaks <- 
    sort(c(-fold_change_thresh, fold_change_thresh,
           seq(0, floor(-x.axis.max/0.2)*0.2, -axisdh),
           seq(0, ceiling(x.axis.max/0.2)*0.2, axisdh)))
  
  if(x.axis.max.nodetect < max(max(xbreaks), abs(min(xbreaks)))){
    x.axis.max.nodetect <- max(max(xbreaks), abs(min(xbreaks))) + 
      (0.5*max(max(xbreaks), abs(min(xbreaks))))
  }
  #xbreaks <- c(xbreaks, x.axis.max.nodetect, -x.axis.max.nodetect)
  xbreaks <- sort(unique(xbreaks))
  
  x.axis <- plot_data$avg_log2FC.x
  x.axis[x.axis == Inf] <- x.axis.max.nodetect
  x.axis[x.axis == -Inf] <- -x.axis.max.nodetect
  x.axis[is.na(x.axis)] <- 0
  
  xlabels <- as.character(xbreaks)
  xlabels[xlabels == x.axis.max.nodetect | xlabels == -x.axis.max.nodetect] <- 
    sprintf("Detected\nonly in\none group")
  
  #mirroring x.axis as y.axis too
  y.axis <- plot_data$avg_log2FC.y
  y.axis[y.axis == Inf] <- x.axis.max.nodetect
  y.axis[y.axis == -Inf] <- -x.axis.max.nodetect
  y.axis[is.na(y.axis)] <- 0
  ybreaks <- xbreaks
  ylabels <- xlabels
  
  #updating plot_data
  plot_data$avg_log2FC.x <- x.axis
  plot_data$avg_log2FC.y <- y.axis
  
  #generating the ggplot object
  scatter_plot <- ggplot(data = plot_data,
                         aes(x = plot_data$avg_log2FC.x, 
                             #x = x.axis, 
                             y = plot_data$avg_log2FC.y, 
                             #y = y.axis, 
                             colour = plot_data$coloured, 
                             text = rownames(plot_data))) + 
    geom_segment(aes(x = -fold_change_thresh, y = fold_change_thresh, 
                     xend = fold_change_thresh, yend = fold_change_thresh), 
                 colour = "grey", size = 0.25, linetype = 2) + 
    geom_segment(aes(x = fold_change_thresh, y = fold_change_thresh, 
                     xend = fold_change_thresh, yend = -fold_change_thresh), 
                 colour = "grey", size = 0.25, linetype = 2) + 
    geom_segment(aes(x = fold_change_thresh, y = -fold_change_thresh, 
                     xend = -fold_change_thresh, yend = -fold_change_thresh), 
                 colour = "grey", size = 0.25, linetype = 2) + 
    geom_segment(aes(x = -fold_change_thresh, y = -fold_change_thresh, 
                     xend = -fold_change_thresh, yend = fold_change_thresh), 
                 colour = "grey", size = 0.25, linetype = 2) + 
    geom_hline(yintercept = 0, colour = "#4d4d4d", size = 0.25) + 
    geom_vline(xintercept = 0, colour = "#4d4d4d", size = 0.25) + 
    geom_abline(slope = 1, intercept = 0, alpha = 0.25)
  #geom_point(shape = 19, size = 1.5) + 
  #adding a regression line and adding pearson correlation values
  if(show_corr){
    #geom_smooth not producing regression line, so calculating coef and adding
    # line manually
    lm_res <- lm(formula = y ~ x, 
                 data = data.frame(x = plot_data$avg_log2FC.x, 
                                   y = plot_data$avg_log2FC.y))
    cor_res <- cor(x = plot_data$avg_log2FC.x, 
                   y = plot_data$avg_log2FC.y, 
                   method = "pearson", use = "complete.obs")
    cor_res <- round(x = cor_res, digits = 3)
    
    scatter_plot <- scatter_plot + 
      geom_abline(slope = coef(lm_res)[['x']], 
                  intercept = coef(lm_res)[['(Intercept)']], 
                  linetype = "dashed", alpha = 0.65) + 
      annotate("text", 
               x = 0+(max(xbreaks)*0.05), 
               #y = max(ybreaks)+0.25,
               y = max(ybreaks)+(0.1*max(ybreaks)), 
               alpha = 0.5, label = paste0("r = ", cor_res),
               hjust = 0, vjust = 1)#, 
    #parse = parse_labels)
    
  }
  scatter_plot <- scatter_plot + 
    geom_point(shape = 19, size = 1.5) + 
    scale_color_manual(
      labels = labels,
      values = values) +
    labs(colour = "Genes", 
         y = paste0("Average log2 fold change", " ", label_y), 
         x = paste0("Average log2 fold change", " ", label_x), 
         title = main_title
    ) + 
    theme_light(base_size = 15) + 
    theme(plot.caption = element_text(hjust = 0, face= "italic"), 
          legend.key.width = unit(0.15,"line")
    )
  
  scatter_plot <- scatter_plot +
    scale_x_continuous(breaks = xbreaks, 
                       #limits = c(min(xbreaks)-1, max(xbreaks)+1), 
                       #limits = c(min(xbreaks)-0.25, max(xbreaks)+0.25), 
                       limits = c(min(xbreaks)+(0.1*min(xbreaks)), 
                                  max(xbreaks)+(0.1*max(xbreaks))), 
                       labels = xlabels) +
    scale_y_continuous(breaks = ybreaks, 
                       #limits = c(min(ybreaks)-1, max(ybreaks)+1), 
                       #limits = c(min(ybreaks)-0.25, max(ybreaks)+0.25), 
                       limits = c(min(ybreaks)+(0.1*min(ybreaks)), 
                                  max(ybreaks)+(0.1*max(ybreaks))), 
                       labels = ylabels) + 
    coord_fixed() + 
    theme(panel.grid.minor.x = element_blank(), 
          #panel.grid.major.x = element_blank(), 
          panel.grid.minor.y = element_blank(), 
          #panel.grid.major.y = element_blank(),
          #axis.line = element_line(colour = 'black')
          panel.grid.major = element_line(colour = "#4d4d4d"),
          axis.ticks = element_line(colour = "#4d4d4d")
    )
  
  #adding quadrant labels
  if(label_quad1 != "" & label_quad3 != ""){
    
    scatter_plot <- scatter_plot + 
      annotate("text", 
               # x = max(xbreaks)+0.25, y = max(ybreaks)+0.25,
               x = max(xbreaks)+(0.1*max(xbreaks)), 
               y = max(ybreaks)+(0.1*max(ybreaks)), 
               alpha = 0.35, label = paste0("Higher in\n", label_quad1), 
               hjust = 1, vjust = 1) + #, 
      #parse = parse_labels) +
      # annotate("text", x = max(xbreaks)+0.25, y = min(ybreaks)-0.25,
      #          alpha = 0.35, label = paste0("Higher in ", label_quad1,
      #                                       ",\nLower in ", label_quad3), 
      #          hjust = 1, vjust = 0) +
      annotate("text", 
               # x = min(xbreaks)-0.25, y = min(ybreaks)-0.25,
               x = min(xbreaks)+(0.1*min(xbreaks)), 
               y = min(ybreaks)+(0.1*min(ybreaks)), 
               alpha = 0.35, label = paste0("Higher in\n", label_quad3), 
               hjust = 0, vjust = 0)#, 
    #parse = parse_labels) #+
    # annotate("text", x = min(xbreaks)-0.25, y = max(ybreaks)+0.25, 
    #          alpha = 0.35, label = paste0("Higher in ", label_quad3,
    #                                       ",\nLower in ", label_quad1), 
    #          hjust = 0, vjust = 1)
  }
  
  scatter_plot <- scatter_plot + ggExtra::removeGrid()
  
  return(scatter_plot)
  
  
}