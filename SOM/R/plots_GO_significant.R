library(readr)
library(plyr)
library(dplyr)
library(RColorBrewer)

##########################Cluster 3

dir <- setwd("/media/mirnis/Seagate\ Backup\ Plus\ Drive/backup/AnastrephaLudens/16.StringTie/SOM/SOM_6,6,rectangular")
samples <- read.table(file.path(dir, "significant_cutree3"), header = FALSE)
files <- file.path(dir, samples$V1)
all(file.exists(files))
GO_cluster3<-lapply(files, read_tsv)
View(GO_cluster3[[1]])
dim(GO_cluster3[[2]])
dim(GO_cluster3[[3]])


ggplot(data=GO_cluster3[[1]], aes(x=GO.ID, y=Significant)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  coord_flip() + scale_fill_brewer(palette="Blues") +
  theme_minimal()

ggplot(data=GO_cluster3[[2]], aes(x=Term, y=Significant)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  coord_flip() + scale_fill_brewer(palette="Blues") +
  theme_minimal()

ggplot(data=GO_cluster3[[3]], aes(x=Term, y=Significant)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  coord_flip() + scale_fill_brewer(palette="Blues") +
  theme_minimal()

##########################Cluster 4

##########################Cluster 5
