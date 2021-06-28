library (gplots)
library (RColorBrewer)
library(viridis)

load("TPM_normalization/headtmap.RData")

coloritos <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(600))
y <- as.matrix(samp2)
hr <- hclust(as.dist(1-cor(t(y), method="pearson")))
hc <- hclust(as.dist(1-cor(y, method="pearson")), method="complete")
mycl  <- cutree(hr, h=max(hr$height)/3); mycolhc <-  topo.colors(length(unique(mycl))); mycolhc <-  mycolhc[as.vector(mycl)]
x<- heatmap.2 (as.matrix(samp2), key=T, symkey=T, trace="none", col=coloritos, Rowv=T, keysize=0.5, cexRow=0.1, cexCol=0.1, density.info="none")
dend.column <- x$colDendrogram

pdf("heatmap_TPM.pdf")

heatmap.2(y,  Rowv=as.dendrogram(hr), Colv=dend.column, col=coloritos,   density.info="none", trace="none", RowSideColors=mycolhc, cexRow=0.1, cexCol=0.7)

dev.off()