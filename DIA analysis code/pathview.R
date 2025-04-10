library(org.Hs.eg.db)  # 如果是小鼠的数据就换成org.Mm.eg.db
library(pathview)


# ----------- 读数据 -------------------
# code for data loading
gene.dat.path <- "./res2_reselect/gene_data_sign.csv"
gene_dat <- read.csv(gene.dat.path)
# ...

gene_dat <- gene_dat[, c("Accession", "FC")]
gene_dat$log2FC <- log2(gene_dat$FC)


gene_dat <- gene_dat[, c("Accession", "log2FC")]

# > head(gene_dat)
# Accession log2FC
#       ...    ...

# ----------- 换名 --------------

gene_map <- AnnotationDbi::select(org.Hs.eg.db, 
                                  keys=gene_dat$Accession,  # TODO: 这里要换成你数据里基因Accession的列
                                  keytype="UNIPROT", 
                                  columns="ENTREZID")

colnames(gene_map)[1] <- "Accession"
entr_dat <- merge(gene_map, gene_dat, by="Accession")
entr_dat <- na.omit(entr_dat)
entr_dat <- entr_dat[!duplicated(entr_dat$ENTREZID), ]

entr_dat <- entr_dat[, c("ENTREZID", "log2FC")]


data = entr_dat$log2FC
names(data) = entr_dat$ENTREZID
head(data)

# ----------- 画图 --------------------
pathview(gene.data = data, 
         limit = list(gene=1, cpd=1),# limit调整颜色bar的上下值
         bins = list(gene = 1000, cpd=1000),
         map.cpdname = T,
         pathway.id = "hsa00020",  # TODO: 换成要画的通路的KEGG ID
         species = "hsa", 
         kegg.native = T,
         split.group = F,
         out.suffix = "citrate cycle_50um",  # TODO：换成匹配要画的通路的suffix
         cex = 0.75, 
         plot.col.key = F,
         new.signature = F
)
