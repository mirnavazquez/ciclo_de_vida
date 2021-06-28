library(DESeq2)
library("devtools")
install_github("kassambara/factoextra")
library("factoextra")

rv <- rowVars(assay(vsd))
select <- order(rv, decreasing = TRUE)[seq_len(min(500, 
                                                   length(rv)))]
pca <- prcomp(t(assay(vsd)[select, ]))
head(pca)
pca$sdev
get_eig(pca)
fviz_eig(pca, addlabels=TRUE, hjust = -0.3) +
  ylim(0, 80)

fviz_eig(pca, choice = "eigenvalue", 
         addlabels=TRUE)
