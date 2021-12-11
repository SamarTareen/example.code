

# function to get cell proportion plots given cluster or other categorisation
# copyright Samar H. K. Tareen

get_cell_counts_props <- function(seurat_object, 
                                  cluster_column = "seurat_clusters", 
                                  category_column = "label", 
                                  #subcategory_column = NA#, 
                                  proportions = T, 
                                  display_stats = T
){
  
  if(!cluster_column %in% colnames(seurat_object@meta.data)){
    stop(paste0("Error: specified cluster column not found in the seurate ", 
                "object. Did you forget to run the respective clustering?"))
  }
  if(!category_column %in% colnames(seurat_object@meta.data)){
    stop(paste0("Error: specified category column not found in the ", 
                "seurate object. Double-check the column names of the meta ", 
                "data."))
  }
  if(length(category_column) > 1){
    stop(paste0("Error: please provide only a single category column."))
  }
  if(length(cluster_column) > 1){
    stop(paste0("Error: please provide only a single cluster column."))
  }
  
  plot_data <- table(seurat_object@meta.data[,category_column], 
                     seurat_object@meta.data[,cluster_column])
  
  #reusing this obsolete part to send by sample proportions to Pierre
  if(proportions){
    #if(proportions & by_sample & !by_cluster){
    plot_data <- (plot_data/rowSums(plot_data))*100
  }
  
  plot_data <- melt(plot_data)
  
  #formatting age_treatment grouping for x-axis
  plot_data$age <- gsub(pattern = "_.*$", 
                        replacement = "", 
                        x = plot_data$Var1)
  plot_data$sample <- gsub(pattern = "^.*_", 
                           replacement = "", 
                           x = plot_data$Var1)
  #plot_data$treatment <- "eGFP"
  plot_data$treatment <- "PHP.GFAP-GFP"
  # plot_data$treatment[
  #   which((plot_data$age %in% c("12m") &
  #            plot_data$sample %in% c("S01", "S03")) |
  #           (!(plot_data$age %in% c("12m")) &
  #              plot_data$sample %in% c("S02", "S04", "S06")) )] <-
  #   "IL2"
  plot_data$treatment[
    which((plot_data$age %in% c("12m") &
             plot_data$sample %in% c("S01", "S03")) |
            (!(plot_data$age %in% c("12m")) &
               plot_data$sample %in% c("S02", "S04", "S06")) )] <-
    "PHP.GFAP-IL2"
  plot_data$age_treatment <- paste0(plot_data$age, 
                                    #"_", 
                                    ".", 
                                    plot_data$treatment)
  
  colplot_list <- list()
  
  #for(i in seq(1, ncol(plot_data))){
  for(i in unique(plot_data$Var2)){
    
    #if(length(dim(plot_data))>2){
    #colplot_data <- melt(plot_data[,i])
    colplot_data <- plot_data[plot_data$Var2==i,]
    
    if(!("Var1" %in% colnames(colplot_data))){
      colplot_data$Var1 <- rownames(colplot_data)
      rownames(colplot_data) <- seq(1, nrow(colplot_data))
    }
    
    col_plot <- 
      ggplot(data = colplot_data, 
             aes(#fill = treatment, 
               #color = treatment, 
               y = value, 
               x = age)) + 
      geom_dotplot(aes(fill = treatment), 
                   binaxis = "y", 
                   #position = "identity", 
                   position = "dodge", 
                   stackdir = "center", 
                   dotsize = 1.5) + 
      stat_summary(aes(fill = treatment), 
                   fun.data=mean_sdl, fun.args = list(mult=1),
                   geom="errorbar", color="black", width=0.2, 
                   position = position_dodge(0.9)) +
      stat_summary(aes(fill = treatment), 
                   fun.y=mean, geom="point", color="black", size = 0.8, 
                   position = position_dodge(0.9)) + 
      # stat_compare_means(method = "t.test", 
      #                    #method = "anova", 
      #                    comparisons = comparisons) + 
      # stat_compare_means(aes(group = age), 
      #                    method = "t.test", 
      #                    comparisons = c("04m", "24m")
      #                    ) + 
      labs(fill = "Treatment", 
           y = "Percentage of cells", 
           x = "Age group", 
           #subtitle = paste0("Cluster: ", as.character(i))
           subtitle = as.character(i)
      ) + 
      expand_limits(y = 0) + 
      theme_light(base_size = 15) + 
      theme(#axis.text.x = element_text(angle = 25, hjust = 1), 
        panel.grid.minor = element_blank(), 
        legend.position = "top")
    
    if(display_stats){
      stat.test.within.age <- colplot_data %>%
        group_by(age) %>%
        t_test(formula = value ~ treatment) %>%
        adjust_pvalue(method = "bonferroni") %>%
        add_significance("p.adj")
      
      stat.test.within.age <- stat.test.within.age %>% 
        add_xy_position(x = "age", dodge = 0.8)
      
      stat.test.within.age$y.position <- stat.test.within.age$y.position + 
        (stat.test.within.age$y.position * 0.05)
      
      col_plot <- col_plot + 
        stat_pvalue_manual(data = stat.test.within.age, 
                           label = "p.adj", tip.length = 0, size = 4.5) #+ 
      #scale_y_continuous(expand = expansion(mult = c(0,0.1)))
      
      stat.test.bw.age <- colplot_data %>%
        t_test(formula = value ~ age) %>%
        adjust_pvalue(method = "bonferroni")
      
      stat.test.bw.age <- stat.test.bw.age %>% 
        add_xy_position(x = "age")
      
      stat.test.bw.age$y.position <- stat.test.bw.age$y.position + 
        (stat.test.bw.age$y.position * 0.2)
      
      col_plot <- col_plot + 
        stat_pvalue_manual(stat.test.bw.age, label = "p.adj", 
                           tip.length = 0.02, step.increase = 0.5, size = 4.5) + 
        scale_y_continuous(expand = expansion(mult = c(0.05,0.1)))
    }
    
    colplot_list[[length(colplot_list)+1]] <- col_plot
    
  }
  
  #paste(toupper(substr(name, 1, 1)), substr(name, 2, nchar(name)), sep="")
  
  return(colplot_list)
  
}