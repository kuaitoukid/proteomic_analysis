library(org.Hs.eg.db)  # Hs for human, Mm for mouse
library(clusterProfiler)
library(tidyverse)
library(enrichplot)

gene.dat <- read.table('./res/1_gene_data_p.txt', header = T)

# >>> filter >>>
# use Group A and B as example
gene.dat <- select(gene.dat, c(Gene, FC_BA, q_AB))
gene.dat <- filter(gene.dat, gene.dat$q_AB < 0.05)
gene.dat$FC_BA <- log2(gene.dat$FC_BA)

# >>> Convert UniprotID to ENTREZID >>>
gene.map <- AnnotationDbi::select(org.Hs.eg.db, 
                                  keys=gene.dat$Gene, 
                                  keytype="UNIPROT", 
                                  columns=("ENTREZID"))
colnames(gene.map)[1] <- "Gene"
entr.dat <- inner_join(gene.map, gene.dat, by="Gene")
entr.dat <- na.omit(entr.dat)

# >>> sort by FoldChange >>>
sort.dat <- entr.dat[order(entr.dat$FC_BA, decreasing = T), ]
gene.lis <- sort.dat[, 3]
names(gene.lis) <- as.character(sort.dat[, 2])

# >>> GSEA >>>
gsea.res <- gseKEGG(gene.lis,
                    organism='hsa',  # http://www.genome.jp/kegg/catalog/org_list.html
                    keyType='kegg',
                    minGSSize=10, 
                    maxGSSize=1000,
                    pvalueCutoff = 0.05)

write.csv(gsea.res, './res/4_GSEA/4_gsea_result.csv', row.names = F)

# >> ridge plot >>
ridge <- ridgeplot(gsea.res, 
                   orderBy = "p.adjust", 
                   decreasing = T)
ggsave("./res/4_GSEA/4_ridgeplot.pdf", ridge, width = 13, height = 13)

# >> gsea plot >>
path <- c("hsa03010", "hsa03040")  # select path from 4_gsea_result.csv
gseaplot <- gseaplot2(gsea.res, path, pvalue_table = F, base_size = 8)
ggsave("./res/4_GSEA/4_gseaplot.pdf", gseaplot, width = 7, height = 5)