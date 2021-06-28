#BiocManager::install("topGO")
library("topGO")
library(ggplot2)

#Universo de genes
geneID2GO1 <- readMappings("data/annotation_03.txt")
geneUniverse1 <- names(geneID2GO1) 

#Clusters
data.val2<-data.val
save(data.val2, file="80,80,dataval2.RData")

enrich_som_cuttree<-function(rango, somData, categoria){
  significant_GO_categoria <- data.frame()
  for (clustNum in c(rango)){  
    message(paste("Running cluster", clustNum))
    sub_cluster <- subset(somData, som_cluster10 == clustNum)
    sub_data <- sub_cluster[,c(2)]
    geneList1 <- factor(as.integer(geneUniverse1 %in% sub_data))
    names(geneList1) <- geneUniverse1
    myGOdata1 <- new("topGOdata", description=paste0("GO Cluster", clustNum), ontology=categoria, allGenes=geneList1,  annot = annFUN.gene2GO, gene2GO = geneID2GO1)
    resultFisher1 <- runTest(myGOdata1, algorithm="weight01", statistic="fisher")
    allGO = usedGO(object = myGOdata1)
    allRes1 <- GenTable(myGOdata1, classicFisher = resultFisher1, orderBy = "resultFisher1", ranksOf = "classicFisher", topNodes = length(allGO))
    #write.table(allRes1, file = paste("Cluster", categoria, "cuttree", clustNum, "all.txt", sep="_"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
    s_GO_categoria<-subset(allRes1, allRes1$classicFisher <= 0.05)
    s_GO_categoria$GO_cluster<-categoria
    s_GO_categoria$cluster<-clustNum
    significant_GO_categoria <- rbind(significant_GO_categoria, s_GO_categoria)
    #write.table(s_GO_categoria, file = paste0("Cluster_", categoria, "_significant_cuttree_", clustNum, ".txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
    showSigOfNodes(myGOdata1, score(resultFisher1), firstSigNodes = 5, useInfo ='all')
  #  printGraph(myGOdata1, resultFisher1, firstSigNodes = 7, fn.prefix = paste0("Cluster_", categoria, clustNum, ".txt"), useInfo = "all", pdfSW = TRUE)
  }
  return(significant_GO_categoria)
}

allClustersBP <- enrich_som_cuttree(1:8, data.val2, "BP")
allClustersCC <- enrich_som_cuttree(1:8, data.val2, "CC")
allClustersMF <- enrich_som_cuttree(1:8, data.val2, "MF")

save(allClustersBP, file="allClustersBP.RData")
save(allClustersCC, file="allClustersCC.RData")
save(allClustersMF, file="allClustersMF.RData")



uno<-subset(allClustersMF, cluster == 1)
dos<-subset(allClustersMF, cluster == 2)
tres<-subset(allClustersMF, cluster == 3)
cuatro<-subset(allClustersMF, cluster == 4)
cinco<-subset(allClustersMF, cluster == 5)
seis<-subset(allClustersMF, cluster == 6)
siete<-subset(allClustersMF, cluster == 7)
ocho<-subset(allClustersMF, cluster == 8)



pdf(file = (paste0("uno_Term_MF.pdf")), width = 8, height = 11)
ggplot(data=uno, aes(x=Term, y=Significant)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  coord_flip()
dev.off()

pdf(file = (paste0("dos_Term_MF.pdf")), width = 8, height = 11)
ggplot(data=dos, aes(x=Term, y=Significant)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  coord_flip()
dev.off()

pdf(file = (paste0("Tres_Term_MF.pdf")), width = 8, height = 11)
ggplot(data=tres, aes(x=Term, y=Significant)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  coord_flip()
dev.off()

pdf(file = (paste0("cuatro_Term_MF.pdf")), width = 8, height = 11)
ggplot(data=cuatro, aes(x=Term, y=Significant)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  coord_flip()
dev.off()

pdf(file = (paste0("cinco_Term_MF.pdf")), width = 8, height = 11)
ggplot(data=cinco, aes(x=Term, y=Significant)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  coord_flip()
dev.off()

pdf(file = (paste0("seis_Term_MF.pdf")), width = 8, height = 11)
ggplot(data=seis, aes(x=Term, y=Significant)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  coord_flip()
dev.off()

pdf(file = (paste0("siete_Term_MF.pdf")), width = 8, height = 11)
ggplot(data=siete, aes(x=Term, y=Significant)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  coord_flip()
dev.off()

pdf(file = (paste0("ocho_Term_MF.pdf")), width = 8, height = 11)
ggplot(data=ocho, aes(x=Term, y=Significant)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  coord_flip()
dev.off()

Extract_Genes_From_GO<-function(rango, somData, enrichout, categoria, expression, IPSannotation){
  GoGene <- data.frame()
  InterProAnnotation <- read.table(IPSannotation, skip=8, stringsAsFactors=FALSE, sep="\t", header=TRUE)
  for (clustNum in c(rango)){  
    message(paste("Running cluster", clustNum))
    sub_cluster <- subset(somData, som_cluster10 == clustNum)
    sub_data <- sub_cluster[,c(2)]
    geneList1 <- factor(as.integer(geneUniverse1 %in% sub_data))
    names(geneList1) <- geneUniverse1
    myGOdata1 <- new("topGOdata", description=paste0("GO Cluster", clustNum), ontology=categoria, allGenes=geneList1,  annot = annFUN.gene2GO, gene2GO = geneID2GO1)
    myterms=enrichout$GO.ID
    mygenes <- genesInTerm(myGOdata1, myterms)
    for (i in 1:length(myterms)) {
      myterm <- myterms[i]
      GoDesc <- unique(enrichout[which(enrichout$GO.ID == myterm),2])
      FullName <- paste0(GoDesc, " (", myterm, ")")
      mygenesforterm <- mygenes[myterm][[1]]
      Gene2Function <- c()
      # Add the PFAM InterProScan annotation. If more than one term is found, they are concatenated
      for (geneGO in mygenesforterm) {
        coordsGeneAnno <- grep(geneGO, InterProAnnotation$Name)
        if (length(coordsGeneAnno) >0) {
          sub <- InterProAnnotation[coordsGeneAnno,]
          PFAM <- paste(geneGO, paste( sub[which(sub$Source ==  "PFAM"), 6], collapse = '|'), sep = '|')
          Gene2Function <- c(Gene2Function, PFAM)
        } else {
          Gene2Function <- c(Gene2Function, geneGO)
        }
      }
      GO <- rep(myterm, length(mygenesforterm))
      GoFullName <- rep(FullName, length(mygenesforterm))
      Cluster <- as.numeric(rep(clustNum, length(mygenesforterm)))
      Gene <- mygenesforterm
      tmp.Go.Gene <- cbind.data.frame(GO, GoFullName, Cluster, Gene, Gene2Function)
      GoGene <- rbind(GoGene, tmp.Go.Gene)
    }
  }
  
  valExp <- expression
  valExp$Gene <- rownames(valExp)
  GoGeneExp <- merge(GoGene, valExp, by="Gene")
  GoGeneExp$Gene <- GoGeneExp$Gene2Function
  GoGeneExp <- GoGeneExp[,-5]
  return(GoGeneExp)
}
# GoGeneExpBP <- Extract_Genes_From_GO(1:8, data.val2, allClustersBP, "BP", samp3, "Anastrepha_ludens_interproscan_annotation.tab")

samp1<-data.val2[,c(2, 9:14)]
samp1<-data.val2[,c(2, 9:14)]
samp3 <- samp1[,-1]
rownames(samp3) <- samp1[,1]

GoGeneExpBP <- Extract_Genes_From_GO(1:8, data.val2, allClustersBP, "BP", samp3)
GoGeneExpCC <- Extract_Genes_From_GO(1:8, data.val2, allClustersCC, "CC", samp3)
GoGeneExpMF <- Extract_Genes_From_GO(1:8, data.val2, allClustersMF, "MF", samp3)

save(GoGeneExpBP, file="BP.RData")
save(GoGeneExpCC, file="CC.RData")
save(GoGeneExpMF, file="MF.RData")



GoGeneExpCC[which(GoGeneExpCC$GO == "GO:0042600" & GoGeneExpCC$Cluster == "8"),]
GoGeneExpMF[which(GoGeneExpMF$GO == "GO:0046332" & GoGeneExpMF$Cluster == "8"),]


unoID<-as.list(unique(uno$GO.ID))
dosID<-as.list(unique(dos$GO.ID))
tresID<-as.list(unique(tres$GO.ID))
cuatroID<-as.list(unique(cuatro$GO.ID))
cincoID<-as.list(unique(cinco$GO.ID))
seisID<-as.list(unique(seis$GO.ID))
sieteID<-as.list(unique(siete$GO.ID))
ochoID<-as.list(unique(ocho$GO.ID))


plotHeatmapGOterms<-function(GOterm){
  message( GOterm)
  Cluster3<-subset(GoGeneExpMF, GO == GOterm)
  GoDesc <- unique(Cluster3$GoFullName)
  gen1<-Cluster3[, -(2:4)]
  colnames(gen1)<-c("Gene", "Adult female", "Adult male", "Egg", "Larvae stage 02", "Larvae stage 03", "Pupae")
  gen4<-c("Gene", "Egg", "Larvae stage 02", "Larvae stage 03", "Pupae", "Adult male", "Adult female")
  gen1 <- gen1[, gen4]
  gen2<-unique(gen1)
  gen3 <- gen2[,-1]
  if(nrow(gen3) > 1 && rowSums(gen3) != 0) {
    rownames(gen3) <- gen2[,1]
    rownames(gen3) <- ifelse(nchar(rownames(gen3)) > 50, paste0(substring(rownames(gen3), 1, 50), "..."), rownames(gen3))
    mat <- t(scale(t(data.matrix(gen3))))
    mat <- mat[rowSums(is.na(mat)) != ncol(mat), ]
    GOtermSan <- gsub(":", "_", GOterm)
    pdf(file = (paste0("Cluster_", GOtermSan, ".pdf")), width = 8, height = 11)
    pheatmap::pheatmap(mat, main= GoDesc, cluster_cols=FALSE, angle_col = "45")
    dev.off()
  } else {
    message("Skipping GO term ", GOterm, ".")
  }
}

lapply(unoID, plotHeatmapGOterms)
lapply(dosID, plotHeatmapGOterms)
lapply(tresID, plotHeatmapGOterms)
lapply(cuatroID, plotHeatmapGOterms)
lapply(cincoID, plotHeatmapGOterms)
lapply(seisID, plotHeatmapGOterms)
lapply(sieteID, plotHeatmapGOterms)
lapply(ochoID, plotHeatmapGOterms)




  save(data.val2, geneID2GO1, geneUniverse1, allClusters, GoGeneExp, samp2, file="DAtos_GO_heatmap.RData")

