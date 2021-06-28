
data.val3<-data.val
names(data.val3)
for (clustNum in c(1:9)){  
  data.val3$cluster[data.val3[,1] == clustNum] <- "subcluster"
  data.val3$cluster[data.val3[,1] != clustNum] <- "other"

  #plot
  
  p <- ggplot(data.val3, aes(PC1, PC2, color = cluster)) 
  p + geom_point(size = I(2), alpha = 0.6) +
    scale_colour_manual(values = c("#cccccc", "#000000")) + 
    theme_bw() + 
    theme(legend.text = element_text(
      size = 30, 
      face = "bold"), 
      text = element_text(size = 30), 
      legend.position = "none")
  ggsave(paste0("Cluster_PCA_cuttree", clustNum, ".pdf"))
}
