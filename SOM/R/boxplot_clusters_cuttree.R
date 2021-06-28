commander<-c("Adult Female", "Adult Male", "Pupa")

for (clustNum in c(1:30)){  
  sub_cluster <- subset(data.val2, som_cluster == clustNum)
  sub_data <- sub_cluster[,c(6:8)] # just the sample types
  colnames(sub_data)<-c("Adult Female", "Adult Male", "Pupa")
  m.data <- melt(sub_data)
  m.data$variable<-factor(m.data$variable, commander)
  p <- ggplot(m.data, aes(x = variable, y =value))
  p + geom_point(alpha = 0.5,position = "jitter", size = 0.5) + 
    geom_boxplot(alpha = 0.75, outlier.size = 0) + 
    theme_bw() + 
    theme(text = element_text(size = 15),
          axis.text.x = element_text(angle = 90, 
                                     vjust = 1)) +
    ggtitle(paste0("Cluster_", clustNum)) +
    xlab("Life stage") +
    ylab("Scaled Gene Expression")
  ggsave(paste0("Cluster_cuttree", clustNum, ".pdf"))
}


