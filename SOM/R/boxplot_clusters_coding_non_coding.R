library(RColorBrewer)

licRNAs <- read.table("data/licRNAs_id_R.txt")
rRNA <- read.table("data/rRNA.txt")
TE <- read.table("data/TE.txt")
tRNA <- read.table("data/tRNA.txt")
data.val4<-data.val


licRNAs$V1<-as.factor(licRNAs$V1)
rRNA$V1<-as.factor(rRNA$V1)
TE$V1<-as.factor(TE$V1)
tRNA$V1<-as.factor(tRNA$V1)

data.val4$gene<-as.factor(data.val4$gene)


for(i in 1:length(data.val4$gene)){
  if (any(as.character(data.val4$gene[i]) == as.character(licRNAs$V1))) {
    data.val4$coding[i] <- c("long_non_coding")
    
  } else if( any(as.character(data.val4$gene[i] == as.character(rRNA$V1)))) {
    data.val4$coding[i] <- c("rRNA")
    
  } else if( any(as.character(data.val4$gene[i] == as.character(TE$V1)))) {
    data.val4$coding[i] <- c("TE")
    
  } else if( any(as.character(data.val4$gene[i] == as.character(tRNA$V1)))) {
    data.val4$coding[i] <- c("tRNA")
    
  } else {
    data.val4$coding[i] <- c("coding")
    
  }
}

head(data.val4)
#commander<-c("Egg", "Larvae stage 2", "Larvae stage 3", "Pupa", "Adult male", "Adult female")
commander<-c("Egg", "Larvae stage 2", "Larvae stage 3", "Pupa")


for (clustNum in c(1:7)){  
  sub_cluster <- subset(data.val4,  som$unit.classif == clustNum & data.val4$coding == "coding" )
  sub_data <- sub_cluster[,c(7:10,16)] # just the sample types
  #colnames(sub_data)<-c("Adult female", "Adult male", "Egg", "Larvae stage 2", "Larvae stage 3", "Pupa", "Sequences")
  colnames(sub_data)<-c("Egg", "Larvae stage 2", "Larvae stage 3", "Pupa", "Sequences")
  m.data <- melt(sub_data)
  m.data$variable<-factor(m.data$variable, commander)
  p <- ggplot(m.data, aes(x = variable, y =value, color=Sequences))
  p + geom_point(alpha = 0.5,position = "jitter", size = 0.5) + 
    geom_boxplot(alpha = 0.75, outlier.size = 0) + 
    theme_bw() + 
    theme(text = element_text(size = 15),
          axis.text.x = element_text(angle = 90, 
                                     vjust = 1)) +
    ggtitle(paste0("CLuster_", clustNum)) +
    xlab("Life stage") +
    ylab("Scaled Gene Expression")
  ggsave(paste0("Cluster_Coding", clustNum, ".pdf"))
}

for (clustNum in c(1:7)){  
  sub_cluster <- subset(data.val4,  som$unit.classif == clustNum & data.val4$coding == "long_non_coding" )
  sub_data <- sub_cluster[,c(7:10,16)] # just the sample types
  #colnames(sub_data)<-c("Adult female", "Adult male", "Egg", "Larvae stage 2", "Larvae stage 3", "Pupa", "Sequences")
  colnames(sub_data)<-c("Egg", "Larvae stage 2", "Larvae stage 3", "Pupa", "Sequences")
  m.data <- melt(sub_data)
  m.data$variable<-factor(m.data$variable, commander)
  p <- ggplot(m.data, aes(x = variable, y =value, color=Sequences))
  p + geom_point(alpha = 0.5,position = "jitter", size = 0.5, color = "#C4961A") + 
    geom_boxplot(alpha = 0.75, outlier.size = 0, fill = "#FFDB6D", color = "#C4961A") + 
    theme_bw() + 
    theme(text = element_text(size = 15),
          axis.text.x = element_text(angle = 90, 
                                     vjust = 1)) +
    ggtitle(paste0("Cluster_", clustNum)) +
    xlab("Life stage") +
    ylab("Scaled Gene Expression") +
  ggsave(paste0("Cluster_long_non_coding", clustNum, ".pdf"))
}
