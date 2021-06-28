library(ggplot2)
library(reshape)
library(kohonen)
library(RColorBrewer)

load("data/SOM1.RData")
mostDEgenes<-allGeneList[c(6, 1, 3)]

which(mostDEgenes$mean==0)
zero<-which(mostDEgenes$mean==0)
mostDEgenes <- mostDEgenes[-c(zero), ]
which(mostDEgenes$mean==0)

sum(is.na(mostDEgenes))

mostDEgene.long <- cast(mostDEgenes, gene ~ stage, value.var = mean, fun.aggregate = "mean")
mostDEgene.long <- as.data.frame(mostDEgene.long)

is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))

mostDEgene.long[is.nan(mostDEgene.long)] <- 0
new_DF <- mostDEgene.long[rowSums(is.na(mostDEgene.long)) > 0,]
new_DF

licRNAs <- read.table("data/licRNAs_id_R.txt")
rRNA <- read.table("data/rRNA.txt")
TE <- read.table("data/TE.txt")
tRNA <- read.table("data/tRNA.txt")

licRNAs$V1<-as.factor(licRNAs$V1)
rRNA$V1<-as.factor(rRNA$V1)
TE$V1<-as.factor(TE$V1)
tRNA$V1<-as.factor(tRNA$V1)

mostDEgene.long$gene<-as.factor(mostDEgene.long$gene)


for(i in 1:length(mostDEgene.long$gene)){
  if (any(as.character(mostDEgene.long$gene[i]) == as.character(licRNAs$V1))) {
    mostDEgene.long$coding[i] <- c("long_non_coding")
    
  } else if( any(as.character(mostDEgene.long$gene[i] == as.character(rRNA$V1)))) {
    mostDEgene.long$coding[i] <- c("rRNA")
    
  } else if( any(as.character(mostDEgene.long$gene[i] == as.character(TE$V1)))) {
    mostDEgene.long$coding[i] <- c("TE")
    
  } else if( any(as.character(mostDEgene.long$gene[i] == as.character(tRNA$V1)))) {
    mostDEgene.long$coding[i] <- c("tRNA")
    
  } else {
    mostDEgene.long$coding[i] <- c("coding")
    
  }
}



mostDEgene.long.cod<-(subset(mostDEgene.long, coding == "coding")) 


scale_data <- as.matrix(t(scale(t(mostDEgene.long.cod[c(2:7)]))))
scale_data[is.nan(scale_data)] <- 0

pca <- prcomp(scale_data, scale=TRUE) 
summary(pca) 


pca.scores <- data.frame(pca$x)
data.val <- cbind(mostDEgene.long.cod, scale_data, pca.scores) 
head(data.val)

colnames(data.val) <- make.unique(names(data.val))

p <- ggplot(data.val, aes(PC1, PC2)) 
p + geom_point()
ggsave("PCA_sin_SOM.pdf")

# Check where the scaled gene expression values are.
names(data.val)
#Subset for SOM in a matrix.  
#som() only works on matrices NOT dataframes
#subset only the scaled gene expression values
som.data <- as.matrix(data.val[,c(9:14)])  

save(som.data, file="SOM_coding.RData")

# Set seed, just make sure you keep the same. 
# Has to do with the randomization process. 
set.seed(1000)

#This is where you change the size of the map
som <- som(som.data, somgrid(150,150,"hexagonal")) 

summary(som)

plot(som, type ="changes")
ggsave("Changes.pdf")

plot(som, type = "codes")
ggsave("Codes.pdf")

plot(som, type = "counts")
ggsave("counts.pdf")

plot(som, type="dist.neighbours")
ggsave("Dist.neighbours.pdf")


## use hierarchical clustering to cluster the codebook vectors
codes<-unlist(som$codes)
codes_som_num<-as.numeric(codes)


som_cluster <- cutree(hclust(dist(codes_som_num)), 8)
# plot these results:
plot(som, type="mapping", bgcol = som_cluster, main = "Clusters") 
add.cluster.boundaries(som, som_cluster) 

som_clusterKey <- data.frame(som_cluster)
som_clusterKey$unit.classif <- c(1:54)

data.val <- cbind(data.val,som$unit.classif,som$distances)
names(data.val)


#Merge data.val with som_clusterKey
##change data.val to match som_cluster key 
names(data.val)[20] <- "unit.classif"

data.val <- merge(data.val, som_clusterKey, by.x = "unit.classif" ) #ignore warning, this is what you want.  You are essentially filling in the empties with the value from som_clusterKey
names(data.val)
dim(data.val)

p <- ggplot(data.val, aes(PC1, PC2, colour=factor(unit.classif))) 
p + geom_point() + theme_bw()
ggsave("Unit.classif.pdf")

p <- ggplot(data.val, aes(PC1, PC2, colour=factor(som_cluster))) 
p + geom_point() + theme_bw()
ggsave("Som_cluster.pdf")

data.val2<-data.val

save(mostDEgenes, mostDEgene.long, scale_data, pca, pca.scores, data.val, som.data, som, codes, codes_som_num, som_cluster, data.val2, file="SOM2_100,100,hexagonal_seed1000.RData")
