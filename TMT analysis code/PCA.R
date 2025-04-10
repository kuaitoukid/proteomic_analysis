library(PCAtools)
library(corrplot)
library(pheatmap)
library(ggforce)
library(tidyverse)

PCA <- function(gene.dat, sample.info, output.dir) {
  # >>> data preparation >>>
  row.names(gene.dat) <- gene.dat$Accession
  gene.dat <- subset(gene.dat, select = -Accession)
  row.names(sample.info) <- sample.info$Sample
  sample.info <- subset(sample.info, select = -Sample)
  
  # >>> correlation >>>
  cor <- cor(gene.dat)
  corplot <- pheatmap(cor)
  ggsave(paste0(output.dir, "corplot.pdf"), corplot, width = 12, height = 6,limitsize = FALSE)
  
  # dendrogram
  dist <- dist(t(gene.dat))
  hc <- hclust(dist)
  pdf(file = paste0(output.dir, "dendrogram.pdf"))
  plot(hc)
  dev.off()
  
  # scree plot
  pca.dat <- pca(gene.dat, metadata = sample.info, scale = TRUE)
  screeplot <- screeplot(pca.dat)
  ggsave(paste0(output.dir, "screeplot.pdf"), screeplot, width = 12, height = 6,limitsize = FALSE)
  
  # PCA
  pca <- biplot(pca.dat, colby = "Group", encircle = T) + 
    # geom_mark_ellipse(aes(fill = sample.info$Group,
    #                       color = sample.info$Group), 
    #                   show.legend=F) +
    theme(legend.position = 'bottom') + 
    coord_equal()
  ggsave(paste0(output.dir, "pca.pdf"), pca, width = 9, height = 6.5)
  
  # loading
  loading <- plotloadings(pca.dat)
  ggsave(paste0(output.dir, "loading.pdf"), loading, width = 12, height = 10)
}