library(tximport)
library(readr)
library(DESeq2)
library("vsn")
library("pheatmap")
library("RColorBrewer")
library(ggplot2)
library(apeglm)
library("genefilter")

dir <- getwd()
samples <- read.table(file.path(dir, "samples.txt"), header = TRUE)
files <- file.path(dir, samples$Name)
all(file.exists(files))
names(files) <-gsub("/t_data.ctab","", samples$Name, perl = T)

tmp <- read_tsv(files[1])
tx2gene <- tmp[, c("t_name", "gene_name")]

txi.hisat <- tximport(files, type = "stringtie", tx2gene = tx2gene)
names(txi.hisat)
head(txi.hisat$counts)

dds <- DESeqDataSetFromTximport(txi.hisat, samples, ~1)
colData(dds)
dds

#Filter genes with very low expression to reduce noise. Here we will remove all genes with lower than 10 reads in total across all samples.
keep = rowSums(counts(dds)) >= 10
dds = dds[keep,]
dds


dds$Etapa <- relevel(dds$Etapa,"egg")
dds$group <- factor(paste0(dds$Etapa))
design(dds) <- ~ group
dds <- DESeq(dds)



pdf(file = "Dispersion.pdf", width = 8, height = 11)
plotDispEsts(dds)
dev.off()

vsd<- varianceStabilizingTransformation(dds, blind = TRUE)
pdf(file = "VSTTransform.pdf", width = 8, height = 11)
meanSdPlot(assay(vsd))
dev.off()
head(assay(vsd), 3)

sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$Etapa, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pdf(file = "HEATMAPsamplesDist.pdf", width = 8, height = 11)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()

pcaData <- plotPCA(vsd, intgroup=c("Etapa"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Etapa)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
ggsave("PCAsamplesGGPLOT.pdf")

topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 20)
mat  <- assay(vsd)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)["Etapa"])
pdf(file = "TopVarSamples.pdf", width = 8, height = 11)
pheatmap(mat, annotation_col = anno)
dev.off()

resultsNames(dds)

dds$group

resEggLs2 <- results(dds, contrast =c("group", "egg", "ls2"), alpha = 0.05 )
resEggLs3 <- results(dds, contrast =c("group", "egg", "ls3"), alpha = 0.05 )
resEggPupa <- results(dds, contrast =c("group", "egg", "pupa"), alpha = 0.05 )
resEggAM <- results(dds, contrast =c("group", "egg",  "AM"), alpha = 0.05 )
resEggAF <- results(dds, contrast =c("group", "egg", "AF"), alpha = 0.05 )

resLs2ls3 <- results(dds, contrast =c("group", "ls2", "ls3"), alpha = 0.05 )
resLs2Pupa <- results(dds, contrast =c("group",  "ls2", "pupa"), alpha = 0.05 )
resLs2AM <- results(dds, contrast =c("group", "ls2", "AM"), alpha = 0.05 )
resLs2AF <- results(dds, contrast =c("group", "ls2", "AF"), alpha = 0.05 )
resLs2Egg <- results(dds, contrast =c("group","ls2", "egg"), alpha = 0.05 )

resLs3ls2 <- results(dds, contrast =c("group", "ls3", "ls2"), alpha = 0.05 )
resLs3Pupa <- results(dds, contrast =c("group", "ls3", "pupa"), alpha = 0.05 )
resLs3AM <- results(dds, contrast =c("group", "ls3", "AM"), alpha = 0.05 )
resLs3AF <- results(dds, contrast =c("group", "ls3", "AF"), alpha = 0.05 )
resLs3Egg <- results(dds, contrast =c("group", "ls3", "egg"), alpha = 0.05 )

resPupals2 <- results(dds, contrast =c("group", "pupa", "ls2"), alpha = 0.05 )
resPupals3 <- results(dds, contrast =c("group", "pupa", "ls3"), alpha = 0.05 )
resPupaAM <- results(dds, contrast =c("group", "pupa", "AM"), alpha = 0.05 )
resPupaAF <- results(dds, contrast =c("group", "pupa", "AF"), alpha = 0.05 )
resPupaEgg <- results(dds, contrast =c("group", "pupa", "egg"), alpha = 0.05 )

resAMls2 <- results(dds, contrast =c("group", "AM", "ls2"), alpha = 0.05 )
resAMls3 <- results(dds, contrast =c("group", "AM", "ls3"), alpha = 0.05 )
resAMAF <- results(dds, contrast =c("group", "AM", "AF"), alpha = 0.05 )
resAMpupa <- results(dds, contrast =c("group", "AM", "pupa"), alpha = 0.05 )
resAMEgg <- results(dds, contrast =c("group", "AM", "egg"), alpha = 0.05 )

resAFls2 <- results(dds, contrast =c("group", "AF", "ls2"), alpha = 0.05 )
resAFls3 <- results(dds, contrast =c("group", "AF",  "ls3"), alpha = 0.05 )
resAFAM <- results(dds, contrast =c("group", "AF", "AM"), alpha = 0.05 )
resAFpupa <- results(dds, contrast =c("group", "AF", "pupa"), alpha = 0.05 )
resAFEgg <- results(dds, contrast =c("group", "AF", "egg"), alpha = 0.05 )

res.list <- list(resEggLs2, resEggLs3, resEggPupa, resEggAM, resEggAF, resLs2ls3, resLs2Pupa, resLs2AM, resLs2AF, resLs2Egg, resLs3ls2, resLs3Pupa, resLs3AM, resLs3AF, resLs3Egg, resPupals2, resPupals3, resPupaAM, resPupaAF, resPupaEgg, resAMls2, resAMls3, resAMAF, resAMpupa, resAMEgg, resAFls2, resAFls3, resAFAM, resAFpupa, resAFEgg)
namesResults <- c("resEggLs2", "resEggLs3", "resEggPupa", "resEggAM", "resEggAF", "resLs2ls3", "resLs2Pupa", "resLs2AM", "resLs2AF", "resLs2Egg", "resLs3ls2", "resLs3Pupa", "resLs3AM", "resLs3AF", "resLs3Egg", "resPupals2", "resPupals3", "resPupaAM", "resPupaAF", "resPupaEgg", "resAMls2", "resAMls3", "resAMAF", "resAMpupa", "resAMEgg", "resAFls2", "resAFls3", "resAFAM", "resAFpupa", "resAFEgg")

res <- list()
resdata <- list()
resdataCut <- list()

for (i in 1:length(res.list)){
  res[[i]]<-res.list[[i]][order(res.list[[i]]$padj),]
  resdata[[i]] <- merge(as.data.frame(res[[i]]), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
  names(resdata[[i]])[1] <- "Gene"
  resdataCut[[i]] <- subset(resdata[[i]][c(1,3,7)], padj <= 0.05 & log2FoldChange < -1 | padj <= 0.05 & log2FoldChange > 1)
  write.table(subset(resdataCut[[i]][c(1,2)]), file= paste( namesResults[i], "Diffexpr-results_DESEQ2.txt", sep = ""), quote = F, row.names = F, col.names = T, sep = "\t")
  
}

### Crear tablas de resultados para suplementarios
for (i in 1:length(res.list)){
  res[[i]]<-res.list[[i]][order(res.list[[i]]$padj),]
  resdata[[i]] <- merge(as.data.frame(res[[i]]), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
  names(resdata[[i]])[1] <- "Gene"
  resdataCut[[i]] <- subset(resdata[[i]][c(1,2,3,4,5,6,7)])
  write.table(subset(resdataCut[[i]]), file= paste( namesResults[i], "Suplementaria1.txt", sep = ""), quote = F, row.names = F, col.names = T, sep = "\t")
  
}


##### Script para hacer la tabla suplementaria
for (i in 1:length(res.list)){
  res[[i]]<-res.list[[i]][order(res.list[[i]]$padj),]
  resdata[[i]] <- merge(as.data.frame(res[[i]]), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
  names(resdata[[i]])[1] <- "Gene"
  resdataCut[[i]] <- subset(resdata[[i]][c(1,3,7)], padj <= 0.05 & log2FoldChange < -1 | padj <= 0.05 & log2FoldChange > 1)
  write.table(subset(resdataCut[[i]]), file= paste( namesResults[i], "Suplementaria.txt", sep = ""), quote = F, row.names = F, col.names = T, sep = "\t")
  
}


save(samples, txi.hisat, dds, vsd, mat, res, resEggLs2, resEggLs3, resEggPupa, resEggAM, resEggAF, resLs2ls3, resLs2Pupa, resLs2AM, resLs2AF, resLs2Egg, resLs3ls2, resLs3Pupa, resLs3AM, resLs3AF, resLs3Egg, resPupals2, resPupals3, resPupaAM, resPupaAF, resPupaEgg, resAMls2, resAMls3, resAMAF, resAMpupa, resAMEgg, resAFls2, resAFls3, resAFAM, resAFpupa, resAFEgg, file="DESeq_Aludens_Genome.RData")






