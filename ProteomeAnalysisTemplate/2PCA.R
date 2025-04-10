library(PCAtools)
library(corrplot)
library(pheatmap)
library(ggforce)
library(org.Hs.eg.db)
library(tidyverse)

output.dir <- "./res/2_PCA/"

# >>> read data >>>
gene.dat <- read.delim("./res/1_gene_data_p.txt")
sample.info <- read.table("./data/2_pca_data_sample_info.txt", sep='\t', header = T, row.names = 1)

# >>> filter >>>
# q < 0.05
gene.dat <- gene.dat[which(gene.dat$q_AB < 0.05 | gene.dat$q_AC < 0.05), ]
gene.dat <- gene.dat[1: 10]

# >>> plot >>>
row.names(gene.dat) <- gene.dat$Gene
gene.dat <- subset(gene.dat, select = -Gene)

# >>> correlation >>>
cor <- cor(gene.dat)
corplot <- pheatmap(cor)
ggsave(paste0(output.dir, "corplot.pdf"), corplot)

# dendrogram
dist <- dist(t(gene.dat))
hc <- hclust(dist)
pdf(file = paste0(output.dir, "dendrogram.pdf"))
plot(hc)
dev.off()

# scree plot
pca.dat <- pca(gene.dat, metadata = sample_info)
screeplot <- screeplot(pca.dat)
ggsave(paste0(output.dir, "screeplot.pdf"), screeplot)

# PCA
pca <- biplot(pca.dat, colby = "Group") + 
  geom_mark_ellipse(aes(fill = sample_info$Group,
                        color = sample_info$Group), 
                    show.legend=F) +
  theme(legend.position = 'bottom') + 
  coord_equal()
ggsave(paste0(output.dir, "pca.pdf"), pca, width = 7, height = 5)

# loading
loading <- plotloadings(pca.dat)
ggsave(paste0(output.dir, "loading.pdf"), loading, width = 7, height = 5)