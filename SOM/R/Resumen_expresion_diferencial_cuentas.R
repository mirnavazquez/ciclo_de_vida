library (gplots)
library (RColorBrewer)
library(viridis)

mostDEgenes<-allGeneList[c(6, 1, 3)]
dim(mostDEgenes)

which(mostDEgenes$mean <=10)
diez<-which(mostDEgenes$mean <=10)
mostDEgenes <- mostDEgenes[-c(diez), ]
which(mostDEgenes$mean <=10)
dim(mostDEgenes)

length(rownames(subset(mostDEgenes, stage == "AM")))
length(rownames(subset(mostDEgenes, stage == "AF")))
length(rownames(subset(mostDEgenes, stage == "egg")))
length(rownames(subset(mostDEgenes, stage == "ls2")))
length(rownames(subset(mostDEgenes, stage == "ls3")))
length(rownames(subset(mostDEgenes, stage == "pupa")))


mostDEgene.long <- cast(mostDEgenes, gene ~ stage, value.var = mean, fun.aggregate = "mean")
mostDEgene.long <- as.data.frame(mostDEgene.long)

is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))

mostDEgene.long[is.nan(mostDEgene.long)] <- 0
head(mostDEgene.long)
as.matrix(mostDEgene.long[c(2:7)])


coloritos <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(600))

heatmap.2  (as.matrix(mostDEgene.long[c(2:7)]), key=T, symkey=T, trace="none",  col=coloritos, dendrogram = c("column"), Rowv=T, keysize=2, cexRow=0.5,  cexCol=0.5)
y <- as.matrix(histericograma)
hr <- hclust(as.dist(1-cor(t(y), method="pearson")))
hc <- hclust(as.dist(1-cor(y, method="pearson")), method="complete")
mycl  <- cutree(hr, h=max(hr$height)/3); mycolhc <-  topo.colors(length(unique(mycl))); mycolhc <-  mycolhc[as.vector(mycl)]
heatmap.2(y,  Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), col=coloritos,  density.info="none", trace="none", RowSideColors=mycolhc,  cexRow=0.1, cexCol=0.7)
heathistericograma  <- heatmap.2 (as.matrix(histericograma), key=T, symkey=T, trace="none", col=coloritos, Rowv=T, keysize=1, cexRow=0.15, cexCol=0.8, density.info="none")
dend.column <- heathistericograma$colDendrogram
heatmap.2(y,  Rowv=as.dendrogram(hr), Colv=dend.column, col=coloritos,   density.info="none", trace="none", RowSideColors=mycolhc, cexRow=0.1, cexCol=0.7)
names=unique(mycl[hr$order])
write.table(names, file="nombres_cluster", quote = F, row.names = F, col.names = F)
write.table(mycl,  file="cluster_unigenes", quote = F, sep = "\t")

names(mycolhc) <- names(mycl); barplot(rep(10, max(mycl)), col=unique(mycolhc[hr$labels[hr$order]]), horiz=F,  names=unique(mycl[hr$order]))
