library(DESeq2)
library("topGO")
library(ggplot2)

#Universo de genes
geneID2GO1 <- readMappings("data/annotation_03.txt")
geneUniverse1 <- names(geneID2GO1) 
###All groups by stage

enrich_all<-function(categoria){
  stages<-levels(vsd$group)
  significant_GO_categoria <- data.frame()
  for (clustNum in stages){  
    colStage <- which(stages %in% clustNum)
    mat <- head( samp2[ order (-samp2[,colStage] ), ], n=500)
    mat_df<-as.data.frame(mat)
    mat_df$gene_symbol <- rownames(mat_df)
    mat_df$gene_symbol
    sub_data <- mat_df$gene_symbol
    geneList1 <- factor(as.integer(geneUniverse1 %in% sub_data))
    names(geneList1) <- geneUniverse1
    myGOdata1 <- new("topGOdata", description=paste0("GO top all variable",clustNum), ontology=categoria, allGenes=geneList1,  annot = annFUN.gene2GO, gene2GO = geneID2GO1)
    resultFisher1 <- runTest(myGOdata1, algorithm="weight01", statistic="fisher")
    allGO = usedGO(object = myGOdata1)
    allRes1 <- GenTable(myGOdata1, classicFisher = resultFisher1, orderBy = "resultFisher1", ranksOf = "classicFisher", topNodes = length(allGO))
    write.table(allRes1, file = paste0("GO_", categoria, "_", clustNum, "_all.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
    s_GO_categoria<-subset(allRes1, allRes1$classicFisher <= 0.05)
    s_GO_categoria$GO_cluster<-categoria
    s_GO_categoria$stage<-clustNum
    significant_GO_categoria <- rbind(significant_GO_categoria, s_GO_categoria)
    write.table(s_GO_categoria, file = paste0("GO_top_all_variable_",categoria,"_significant_", clustNum, ".txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
    showSigOfNodes(myGOdata1, score(resultFisher1), firstSigNodes = 5, useInfo ='all')
    printGraph(myGOdata1, resultFisher1, firstSigNodes = 7, fn.prefix = paste0("Cluster_", categoria, clustNum, ".txt"), useInfo = "all", pdfSW = TRUE)
  }
  
  return(significant_GO_categoria)
}



enrich_CC<-enrich_all("CC")

Egg_CC<-subset(enrich_CC, stage == "egg")
AM_CC<-subset(enrich_CC, stage == "AM")
AF_CC<-subset(enrich_CC, stage == "AF")
Pupa_CC<-subset(enrich_CC, stage == "pupa")
Ls2_CC<-subset(enrich_CC, stage == "ls2")
Ls3_CC<-subset(enrich_CC, stage == "ls3")


plot_box<-function(subData, nombre){
  ggplot(data=subData, aes(x=Term, y=Significant)) + 
    geom_bar(stat="identity", position=position_dodge()) +
    coord_flip()
  ggsave(paste0(nombre, ".pdf"))
}

plot_box(Egg_CC, "EggCC")
plot_box(AM_CC, "AM_CC")
plot_box(AF_CC, "AF_CC")
plot_box(Pupa_CC, "Pupa_CC")
plot_box(Ls2_CC, "Ls2_CC")
plot_box(Ls3_CC, "Ls3_CC")


Extract_Genes_From_GO<-function(enrichout, categoria, expression){
  stages<-levels(vsd$group)
  GoGene <- data.frame()
  for (clustNum in stages){  
    colStage <- which(stages %in% clustNum)
    mat <- head( samp2[ order (-samp2[,colStage] ), ], n=500)
    mat_df<-as.data.frame(mat)
    mat_df$gene_symbol <- rownames(mat_df)
    mat_df$gene_symbol
    sub_data <- mat_df$gene_symbol
    geneList1 <- factor(as.integer(geneUniverse1 %in% sub_data))
    names(geneList1) <- geneUniverse1
    myGOdata1 <- new("topGOdata", description=paste0("GO Cluster", clustNum), ontology=categoria, allGenes=geneList1,  annot = annFUN.gene2GO, gene2GO = geneID2GO1)
    myterms=enrichout$GO.ID
    mygenes <- genesInTerm(myGOdata1, myterms)
    for (i in 1:length(myterms))
    {
      myterm <- myterms[i]
      mygenesforterm <- mygenes[myterm][[1]]
      GO <- rep(myterm, length(mygenesforterm))
      Gene <- mygenesforterm
      tmp.Go.Gene <- cbind.data.frame(GO,  Gene)
      GoGene <- rbind(GoGene, tmp.Go.Gene)
    }
  }
  
  valExp <- expression
  valExp$Gene <- rownames(valExp)
  GoGeneExp <- merge(GoGene, valExp, by="Gene")
  return(GoGeneExp)
}

GoGeneExp <- Extract_Genes_From_GO(enrich_CC, "CC", samp2)


dim(GoGeneExp)

Egg_CCID<-as.list(unique(Egg_CC$GO.ID))
AM_CCID<-as.list(unique(AM_CC$GO.ID))
AF_CCID<-as.list(unique(AF_CC$GO.ID))
Pupa_CCID<-as.list(unique(Pupa_CC$GO.ID))
Ls2_CCID<-as.list(unique(Ls2_CC$GO.ID))
Ls3_CCID<-as.list(unique(Ls3_CC$GO.ID))



plotHeatmapGOterms<-function(GOterm){
  Cluster3<-subset(GoGeneExp, GO == GOterm)
  message( GOterm)
  gen1<-Cluster3[, -(2)]
  gen2<-unique(gen1)
  gen3 <- gen2[,-1]
  if(nrow(gen3) > 1) {
    rownames(gen3) <- gen2[,1]
    mat <- t(scale(t(data.matrix(gen3))))
    mat <- mat[rowSums(is.na(mat)) != ncol(mat), ]
    GOtermSan <- gsub(":", "_", GOterm)
    pdf(file = (paste0("Cluster_", GOtermSan, ".pdf")), width = 8, height = 11)
    pheatmap::pheatmap(mat, main= GOterm)
    dev.off()
  } else {
    message("Skipping GO term ", GOterm, ".")
  }
}



lapply(Egg_CCID, plotHeatmapGOterms)
lapply(AM_CCID, plotHeatmapGOterms)
lapply(AF_CCID, plotHeatmapGOterms)
lapply(Pupa_CCID, plotHeatmapGOterms)
lapply(Ls2_CCID, plotHeatmapGOterms)
lapply(Ls3_CCID, plotHeatmapGOterms)



save(GoGeneExp, enrich_CC, Egg_CCID, Egg_CC, file ="all_GO_TMP.RData")


