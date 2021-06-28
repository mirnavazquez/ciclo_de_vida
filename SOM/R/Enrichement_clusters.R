#BiocManager::install("topGO")
library("topGO")

#Universo de genes
geneID2GO1 <- readMappings("data/annotation_03.txt")
geneUniverse1 <- names(geneID2GO1) 

#Clusters
data.val2<-data.val
head(data.val2)

for (clustNum in c(1:9)){  
  sub_cluster <- subset(data.val2, som$unit.classif == clustNum)
  sub_data <- sub_cluster[,c(2)]
  geneList1 <- factor(as.integer(geneUniverse1 %in% sub_data))
  names(geneList1) <- geneUniverse1
  myGOdata1 <- new("topGOdata", description=paste0("GO Cluster", clustNum), ontology="BP", allGenes=geneList1,  annot = annFUN.gene2GO, gene2GO = geneID2GO1)
  resultFisher1 <- runTest(myGOdata1, algorithm="weight01", statistic="fisher")
  allGO = usedGO(object = myGOdata1) 
  allRes1 <- GenTable(myGOdata1, classicFisher = resultFisher1, orderBy = "resultFisher1", ranksOf = "classicFisher", topNodes = length(allGO))
  write.table(allRes1, file = paste0("Cluster_BP", clustNum, "_all.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  significant_GO<-subset(allRes1, allRes1$classicFisher <= 0.05)
  write.table(significant_GO, file = paste0("Cluster_BP_significant", clustNum, ".txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
}

for (clustNum in c(1:9)){  
  sub_cluster <- subset(data.val2, som$unit.classif == clustNum)
  sub_data <- sub_cluster[,c(2)]
  geneList1 <- factor(as.integer(geneUniverse1 %in% sub_data))
  names(geneList1) <- geneUniverse1
  myGOdata1 <- new("topGOdata", description=paste0("GO Cluster", clustNum), ontology="CC", allGenes=geneList1,  annot = annFUN.gene2GO, gene2GO = geneID2GO1)
  resultFisher1 <- runTest(myGOdata1, algorithm="weight01", statistic="fisher")
  allGO = usedGO(object = myGOdata1) 
  allRes1 <- GenTable(myGOdata1, classicFisher = resultFisher1, orderBy = "resultFisher1", ranksOf = "classicFisher", topNodes = length(allGO))
  write.table(allRes1, file = paste0("Cluster_CC", clustNum, "_all.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  significant_GO<-subset(allRes1, allRes1$classicFisher <= 0.05)
  write.table(significant_GO, file = paste0("Cluster_CC_significant", clustNum, ".txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
}

for (clustNum in c(1:9)){  
  sub_cluster <- subset(data.val2, som$unit.classif == clustNum)
  sub_data <- sub_cluster[,c(2)]
  geneList1 <- factor(as.integer(geneUniverse1 %in% sub_data))
  names(geneList1) <- geneUniverse1
  myGOdata1 <- new("topGOdata", description=paste0("GO Cluster", clustNum), ontology="MF", allGenes=geneList1,  annot = annFUN.gene2GO, gene2GO = geneID2GO1)
  resultFisher1 <- runTest(myGOdata1, algorithm="weight01", statistic="fisher")
  allGO = usedGO(object = myGOdata1) 
  allRes1 <- GenTable(myGOdata1, classicFisher = resultFisher1, orderBy = "resultFisher1", ranksOf = "classicFisher", topNodes = length(allGO))
  write.table(allRes1, file = paste0("Cluster_MF", clustNum, "_all.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  significant_GO<-subset(allRes1, allRes1$classicFisher <= 0.05)
  write.table(significant_GO, file = paste0("Cluster_MF_significant", clustNum, ".txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
}
