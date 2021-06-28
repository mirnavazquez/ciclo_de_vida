library(readr)
library(plyr)
library(dplyr)

dir <- getwd()
samples <- read.table(file.path(dir, "samples2.txt"), header = TRUE)
files <- file.path(dir, samples$Name)
all(file.exists(files))

##https://stackoverflow.com/questions/11433432/how-to-import-multiple-csv-files-at-once
normalized_counts = lapply(files, read_tsv) %>% bind_rows()

head(normalized_counts)

normalized_counts$stage <- ifelse(grepl("05", normalized_counts$Name, ignore.case = T), "egg", 
                                  ifelse(grepl("11", normalized_counts$Name, ignore.case = T), "egg", 
                                         ifelse(grepl("17", normalized_counts$Name, ignore.case = T), "egg",
                                                ifelse(grepl("06", normalized_counts$Name, ignore.case = T), "ls2",
                                                       ifelse(grepl("12", normalized_counts$Name, ignore.case = T), "ls2",
                                                              ifelse(grepl("18", normalized_counts$Name, ignore.case = T), "ls2",
                                                                     ifelse(grepl("07", normalized_counts$Name, ignore.case = T), "ls3",
                                                                            ifelse(grepl("19", normalized_counts$Name, ignore.case = T), "ls3",
                                                                                   ifelse(grepl("08", normalized_counts$Name, ignore.case = T), "pupa",
                                                                                          ifelse(grepl("14", normalized_counts$Name, ignore.case = T), "pupa",
                                                                                                 ifelse(grepl("20", normalized_counts$Name, ignore.case = T), "pupa",
                                                                                                        ifelse(grepl("09", normalized_counts$Name, ignore.case = T), "AM",
                                                                                                               ifelse(grepl("15", normalized_counts$Name, ignore.case = T), "AM",
                                                                                                                      ifelse(grepl("21", normalized_counts$Name, ignore.case = T), "AM",
                                                                                                                             ifelse(grepl("10", normalized_counts$Name, ignore.case = T), "AF",
                                                                                                                                    ifelse(grepl("16", normalized_counts$Name, ignore.case = T), "AF",
                                                                                                                                           ifelse(grepl("22", normalized_counts$Name, ignore.case = T), "AF","unknown")))))))))))))))))

tail(normalized_counts)

allGenes <- rbind(resEggLs2, resEggLs3, resEggPupa, resEggAM, resEggAF, resLs2ls3, resLs2Pupa, resLs2AM, resLs2AF, resLs2Egg, resLs3ls2, resLs3Pupa, resLs3AM, resLs3AF, resLs3Egg, resPupals2, resPupals3, resPupaAM, resPupaAF, resPupaEgg, resAMls2, resAMls3, resAMAF, resAMpupa, resAMEgg, resAFls2, resAFls3, resAFAM, resAFpupa, resAFEgg)
head(allGenes)
dim(allGenes)

allGenesITAG<-rownames(allGenes)
#allGenesITAG <- allGenes[,1]
length(allGenesITAG)
#Remove duplicates
allGenesITAG <- as.factor(unique(allGenesITAG))
length(allGenesITAG)
class(allGenesITAG)

#make an empty table to hold all the genes
allGeneList <- data.frame(t(rep(NA,7)))
colnames(allGeneList) <- c("type", "genotype", "N", "mean", "sd", "se", "gene")
allGeneList <- allGeneList[-1,] #remove first row

# Takes some time to run 

for(GENE in allGenesITAG) {
  
  if(length(grep(GENE, normalized_counts$Gene_ID)) < 1){ #this is just making sure that the list of sig genes
    next; 
  }
  
  geneData <- subset(normalized_counts, grepl(GENE, normalized_counts$Gene_ID))
  
  sumGraph <- ddply(geneData, "stage", summarise,
                    N    = length(TPM),
                    mean = mean(TPM),  
                    sd   = sd(TPM),
                    se   = sd / sqrt(N))
  
  sumGraph$gene <- GENE
  
  allGeneList  <- rbind(allGeneList, sumGraph) #bind together all the new rows per loop. 
}

dim(allGeneList)
head(allGeneList)

save(normalized_counts, allGenes, allGeneList, file="ExpressionTPM.RData")

mostDEgenes<-allGeneList[c(6, 1, 3)]

which(mostDEgenes$mean <=10)
diez<-which(mostDEgenes$mean <=10)
mostDEgenes <- mostDEgenes[-c(diez), ]
which(mostDEgenes$mean <=10)
dim(mostDEgenes)

length(mostDEgenes$gene)
length(unique(mostDEgenes$gene))


subAM<-subset(mostDEgenes, stage == "AM" )
length(unique(subAM$gene))

subegg<-subset(mostDEgenes, stage == "egg" )
length(unique(subegg$gene))

subls2<-subset(mostDEgenes, stage == "ls2" )
length(unique(subls2$gene))

subls3<-subset(mostDEgenes, stage == "ls3" )
length(unique(subls3$gene))

subpupa<-subset(mostDEgenes, stage == "pupa" )
length(unique(subpupa$gene))

subAF<-subset(mostDEgenes, stage == "AF" )
length(unique(subAF$gene))


library(reshape)

mostDEgene.long <- cast(mostDEgenes, gene ~ stage, value.var = mean, fun.aggregate = "mean")
mostDEgene.long <- as.data.frame(mostDEgene.long)

is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))

mostDEgene.long[is.nan(mostDEgene.long)] <- 0
class(mostDEgene.long)

samp2 <- mostDEgene.long[,-1]
rownames(samp2) <- mostDEgene.long[,1]

dim(samp2)

save(samp2, file="headtmap.RData")

library(pheatmap)

colnames(samp2)<-c("Adult female", "Adult male", "Egg", "Larvae stage 02", "Larvae stage 03", "Pupae")
head(samp2)
samp1<-c("Egg", "Larvae stage 02", "Larvae stage 03", "Pupae", "Adult male", "Adult female")
samp1 <- samp2[, samp1]
head(samp1)

pheatmap::pheatmap(as.matrix(samp1), scale = "row", cluster_cols=FALSE,  angle_col = c("0"))
