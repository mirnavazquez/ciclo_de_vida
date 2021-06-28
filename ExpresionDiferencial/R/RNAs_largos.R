library(RColorBrewer)

licRNAs <- read.table("../SOM/data/licRNAs_id_R.txt")
rRNA <- read.table("../SOM/data/rRNA.txt")
TE <- read.table("../SOM/data/TE.txt")
tRNA <- read.table("../SOM/data/tRNA.txt")



for(i in 1:length(rownames(samp2))){
  if (any(as.character(rownames(samp2)[i]) == as.character(licRNAs$V1))) {
    samp2$coding[i] <- c("long_non_coding")
    
  } else if( any(as.character(rownames(samp2)[i] == as.character(rRNA$V1)))) {
    samp2$coding[i] <- c("rRNA")
    
  } else if( any(as.character(rownames(samp2)[i] == as.character(TE$V1)))) {
    samp2$coding[i] <- c("TE")
    
  } else if( any(as.character(rownames(samp2)[i] == as.character(tRNA$V1)))) {
    samp2$coding[i] <- c("tRNA")
    
  } else {
    samp2$coding[i] <- c("coding")
    
  }
}


colnames(samp2)<-c("Adult female", "Adult male", "Egg", "Larvae stage 02", "Larvae stage 03", "Pupae", "coding")
head(samp2)
samp1<-c("Egg", "Larvae stage 02", "Larvae stage 03", "Pupae", "Adult male", "Adult female", "coding")
samp1 <- samp2[, samp1]

subset_long_non_coding<-subset(samp1, coding == "long_non_coding")
dim(subset_long_non_coding)
subset_rRNA<-subset(samp1, coding == "rRNA")
dim(subset_rRNA)
subset_TE<-subset(samp1, coding == "TE")
dim(subset_TE)
subset_tRNA<-subset(samp1, coding == "tRNA")
dim(subset_tRNA)

collnrna <- colorRampPalette(brewer.pal(10, "BrBG"))(256)
coLTE <- colorRampPalette(brewer.pal(10, "PuOr"))(256)

pheatmap::pheatmap(as.matrix(subset_long_non_coding[,-7]), scale = "row", col=collnrna, cluster_cols=FALSE,  angle_col = c("0"))
pheatmap::pheatmap(as.matrix(subset_TE[,-7]), scale = "row", col=coLTE, cluster_cols=FALSE, angle_col = c("0"))



