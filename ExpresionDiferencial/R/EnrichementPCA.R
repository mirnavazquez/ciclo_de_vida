library(DESeq2)
library("topGO")
library(ggplot2)

#Universo de genes
geneID2GO1 <- readMappings("data/annotation_03.txt")
geneUniverse1 <- names(geneID2GO1) 


####All groups

enrich500<-function(categoria){
  topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 500)
  mat  <- assay(vsd)[ topVarGenes, ]
  mat_df<-as.data.frame(mat)
  mat_df$gene_symbol <- rownames(mat_df)
  sub_data <- mat_df$gene_symbol
  geneList1 <- factor(as.integer(geneUniverse1 %in% sub_data))
  names(geneList1) <- geneUniverse1
  myGOdata1 <- new("topGOdata", description=paste0("GO top 500 variable", categoria), ontology=categoria, allGenes=geneList1,  annot = annFUN.gene2GO, gene2GO = geneID2GO1)
  resultFisher1 <- runTest(myGOdata1, algorithm="weight01", statistic="fisher")
  allGO = usedGO(object = myGOdata1) 
  allRes1 <- GenTable(myGOdata1, classicFisher = resultFisher1, orderBy = "resultFisher1", ranksOf = "classicFisher", topNodes = length(allGO))
  write.table(allRes1, file =paste0( "GO_top_500_variable_",categoria,".txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  significant_GO_categoria<-subset(allRes1, allRes1$classicFisher <= 0.05)
  significant_GO_categoria$GO_cluster<-categoria
  write.table(significant_GO_categoria, file = paste0("GO_top_500_variable_",categoria,"_significant.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE) 
  return(significant_GO_categoria)
}


ggplot(data=BP, aes(x=Term, y=Significant)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  coord_flip()  

ggplot(data=MF, aes(x=Term, y=Significant)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  coord_flip() 

ggplot(data=CC, aes(x=Term, y=Significant)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  coord_flip() 
