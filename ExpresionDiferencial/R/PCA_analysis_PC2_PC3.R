library(ggplot2)
library(ggrepel)
library(DESeq2)

plotPCA.san <- function (object, intgroup=c("Etapa"),  ntop = 500, returnData = FALSE) 
{
  rv <- rowVars(assay(vsd))
  select <- order(rv, decreasing = TRUE)[seq_len(min(500, 
                                                     length(rv)))]
  pca <- prcomp(t(assay(vsd)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all("Etapa" %in% names(colData(vsd)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(vsd)[, "Etapa", drop = FALSE])
  if (length("Etapa") > 1) {
    group <- factor(apply(intgroup.df, 1, paste, collapse = " : "))
  } else {
    group <- colData(vsd)[["Etapa"]]
  }
  
  ## Select the PCAs and percentVar that you like instead of 1 and 2
  d <- data.frame(PC2 = pca$x[, 2], PC3 = pca$x[, 3], group = group, 
                  intgroup.df, name = colData(vsd)[,1])
  if (FALSE) {
    attr(d, "percentVar") <- percentVar[2:3]
    return(d)
  }
  pdf(file = "PC2vsPC3.pdf")
  ggplot(data = d, aes_string(x = "PC2", y = "PC3", color = "group")) + 
    geom_point(size = 3) + 
    xlab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) + 
    ylab(paste0("PC3: ", round(percentVar[3] * 100), "% variance")) + 
    coord_fixed()  
  dev.off()
}

