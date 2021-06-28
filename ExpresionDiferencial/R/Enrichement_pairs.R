library(DESeq2)
library("topGO")
library(ggplot2)
library(dplyr)


resEggLs2_up <- as.data.frame(subset(resEggLs2[c(1,2,6)], padj <= 0.05 & log2FoldChange > 1))
dim(resEggLs2_up)
resEggLs2_down <- as.data.frame(subset(resEggLs2[c(1,2,6)], padj <= 0.05 & log2FoldChange < -1 ))
dim(resEggLs2_down)
resLs2Egg_up <- as.data.frame(subset(resLs2ls3[c(1,2,6)], padj <= 0.05 & log2FoldChange > 1))
dim(resLs2Egg_up)
resLs2Egg_down <- as.data.frame(subset(resLs2ls3[c(1,2,6)], padj <= 0.05 & log2FoldChange < -1 ))
dim(resLs2Egg_down)
resLs3Pupa_up <- as.data.frame(subset(resLs3Pupa[c(1,2,6)], padj <= 0.05 & log2FoldChange > 1))
dim(resLs3Pupa_up)
resLs3Pupa_down <- as.data.frame(subset(resLs3Pupa[c(1,2,6)], padj <= 0.05 & log2FoldChange < -1))
dim(resLs3Pupa_down)
resPupaAM_up <- as.data.frame(subset(resPupaAM[c(1,2,6)], padj <= 0.05 & log2FoldChange > 1))
dim(resPupaAM_up)
resPupaAM_down <- as.data.frame(subset(resPupaAM[c(1,2,6)], padj <= 0.05 & log2FoldChange < -1))
dim(resPupaAM_down)
resPupaAF_up <- as.data.frame(subset(resPupaAF[c(1,2,6)],  padj <= 0.05 & log2FoldChange > 1))
dim(resPupaAF_up)
resPupaAF_down <- as.data.frame(subset(resPupaAF[c(1,2,6)], padj <= 0.05 & log2FoldChange < -1 ))
dim(resPupaAF_down)
resAFAM_up <- as.data.frame(subset(resAFAM[c(1,2,6)],  padj <= 0.05 & log2FoldChange > 1))
dim(resAFAM_up)
resAFAM_down <- as.data.frame(subset(resAFAM[c(1,2,6)], padj <= 0.05 & log2FoldChange < -1 ))
dim(resAFAM_down)




#Universo de genes
geneID2GO1 <- readMappings("data/annotation_03.txt")
geneUniverse1 <- names(geneID2GO1) 
###All groups by stage 

enrich_all<-function(lista, nombre, categoria){
    significant_GO_categoria <- data.frame()
    sub_data <- rownames(lista)
    geneList1 <- factor(as.integer(geneUniverse1 %in% sub_data))
    names(geneList1) <- geneUniverse1
    myGOdata1 <- new("topGOdata", description=paste0("GO top all variable",nombre), ontology=categoria, allGenes=geneList1,  annot = annFUN.gene2GO, gene2GO = geneID2GO1)
    resultFisher1 <- runTest(myGOdata1, algorithm="weight01", statistic="fisher")
    allGO = usedGO(object = myGOdata1)
    allRes1 <- GenTable(myGOdata1, classicFisher = resultFisher1, orderBy = "resultFisher1", ranksOf = "classicFisher", topNodes = length(allGO))
    write.table(allRes1, file = paste0("GO_", categoria, "_", nombre, "_all.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
    s_GO_categoria<-subset(allRes1, allRes1$classicFisher <= 0.01)
    s_GO_categoria$GO_cluster<-categoria
    s_GO_categoria$UpDown<-nombre
    significant_GO_categoria <- rbind(significant_GO_categoria, s_GO_categoria)
    write.table(s_GO_categoria, file = paste0("GO_top_all_variable_",categoria,"_significant_", nombre, ".txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
    showSigOfNodes(myGOdata1, score(resultFisher1), firstSigNodes = 5, useInfo ='all')
    printGraph(myGOdata1, resultFisher1, firstSigNodes = 7, fn.prefix = paste0("Cluster_", categoria, nombre, ".txt"), useInfo = "all", pdfSW = TRUE)
  return(significant_GO_categoria)
}

resAFAM_up_BP_List<-enrich_all(resAFAM_up, "AFAMUP", "BP")
resAFAM_up_CC_List<-enrich_all(resAFAM_up, "AFAMUP", "CC")
resAFAM_up_MF_List<-enrich_all(resAFAM_up, "AFAMUP", "MF")

resAFAM_down_BP_List<-enrich_all(resAFAM_down, "AFAMDown", "BP")
resAFAM_down_CC_List<-enrich_all(resAFAM_down, "AFAMDown", "CC")
resAFAM_down_MF_List<-enrich_all(resAFAM_down, "AFAMDown", "MF")


AFAMBP<-rbind(resAFAM_up_BP_List, resAFAM_down_BP_List)
AFAMCC<-rbind(resAFAM_up_CC_List, resAFAM_down_CC_List)
AFAMMF<-rbind(resAFAM_up_MF_List, resAFAM_down_MF_List)


plot_box<-function(subData, nombre){
  ggplot(data=subData, aes(x=Term, y=log(subData$Significant), fill=UpDown)) + 
    geom_bar(stat="identity", position=position_dodge()) +
    coord_flip() + scale_fill_brewer(palette="Paired")
  ggsave(paste0(nombre, ".pdf"))
}

plot_box(AFAMBP, "AFAMBP")
plot_box(AFAMCC, "AFAMCC")
plot_box(AFAMMF, "AFAMMF")







