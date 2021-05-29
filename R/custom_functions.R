library(reshape2)
library(ggplot2)

##########################
# User defined functions #
##########################

fix_col_names <- function(subacc) {
  col_names <- colnames(subacc)
  # Remove the prefix with the experiment name in all the columns
  new_col_names <- sub("RNA_Seq_v2_mRNA_median_all_sample_Zscores_", "", col_names)
  # Set the new column names
  colnames(subacc) <- new_col_names
  return(subacc)
}

reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
  return (cormat)
}

create_ggheatmap <- function(metled_corr, title) {
  ggheatmap <- ggplot(metled_corr, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name=paste(title,"Pearson","Correlation",sep="\n")) +
    theme_minimal()+ # minimal theme
    coord_fixed() +
    geom_text(aes(Var2, Var1, label = sprintf("%0.2f", round(value, digits = 2))), 
              color = "black", 
              size = 3) +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 1, 
                                 size = 8, hjust = 1),
      axis.text.y = element_text(size = 8),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      legend.justification = c(1, 0),
      legend.position = c(0.6, 0.7),
      legend.direction = "horizontal")+
    guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                                 title.position = "top", title.hjust = 0.5))
  return (ggheatmap)
}

calculate_corr_and_heatmap <- function(subacc, title, gene_signature) {
  # Find genes in common
  genes_in_common <- intersect(gene_signature,colnames(subacc))
  # length(genes_in_common)
  cols_subset <- c(genes_in_common, c("OS_MONTHS","DSS_MONTHS","PFS_MONTHS"))
  # length(cols_subset)
  
  mycors <- cor(as.matrix(subacc[,cols_subset]))
  # Reorder the correlation matrix
  cormat <- reorder_cormat(mycors)
  # cormat <- mycors
  cormat[lower.tri(cormat)] <- NA
  # Melt the correlation matrix
  melted_cormat <- melt(cormat, na.rm = TRUE)
  ggheatmap <- create_ggheatmap(melted_cormat,title)
  return (ggheatmap)
}
