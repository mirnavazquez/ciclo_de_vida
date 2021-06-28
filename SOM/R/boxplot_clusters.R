commander<-c("Egg", "Larvae stage 2", "Larvae stage 3", "Pupa")

for (clustNum in c(1:100)){  
  sub_cluster <- subset(data.val2, som$unit.classif == clustNum)
  sub_data <- sub_cluster[,c(7:10)] # just the sample types
  colnames(sub_data)<-c("Egg", "Larvae stage 2", "Larvae stage 3", "Pupa")
  m.data <- melt(sub_data)
  m.data$variable<-factor(m.data$variable, commander)
  p <- ggplot(m.data, aes(x = variable, y =value))
  p + geom_point(alpha = 0.5,position = "jitter", size = 0.5) + 
    geom_boxplot(alpha = 0.75, outlier.size = 0) + 
    theme_bw() + 
    theme(text = element_text(size = 15),
          axis.text.x = element_text(angle = 90, 
                                     vjust = 1)) +
    ggtitle(paste0("CLuster_", clustNum)) +
    xlab("Life stage") +
    ylab("Scaled Gene Expression")
  ggsave(paste0("Cluster_", clustNum, ".pdf"))
}
