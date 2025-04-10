library(ggplot2)  # 柱状图和点状图
library(stringr)  # 基因ID转换
library(enrichplot)  # GO,KEGG,GSEA
library(clusterProfiler)  # GO,KEGG,GSEA
library(GOplot)  # 弦图，弦表图，系统聚类图
library(DOSE)
library(ggnewscale)
library(topGO)  # 绘制通路网络图
library(circlize)  # 绘制富集分析圈图
library(ComplexHeatmap)  # 绘制图例
library(org.Hs.eg.db)
library(tidyverse)
library(ggforce)

GO_database <- org.Hs.eg.db
KEGG_database <- "hsa"

output.dir <- "./res/5_GO/"
output.fileformat <- c(".png", ".pdf")

# >>> read file >>>
gene.dat <- read.csv("./res/1_gene_data_p.txt", header = T)

# >>> filter >>>
# use Group A and B as example
gene.dat <- select(gene.dat, c(Gene, FC_BA, q_AB))
gene.dat <- filter(gene.dat, gene.dat$q_AB < 0.05 & (gene.dat$FC_BA > 1.3 | gene.dat$FC_BA < 0.8))
colnames(gene.dat) <- c("Gene", "FoldChange", "qvalue")

# >>> Convert UniprotID to ENTREZID >>>
gene.map <- AnnotationDbi::select(GO_database, 
                                  keys=gene.dat$Gene, 
                                  keytype="UNIPROT", 
                                  columns=c("ENTREZID", "SYMBOL"))

colnames(gene.map)[1] <- "Gene"
entr.dat <- inner_join(gene.map, gene.dat, by="Gene")
entr.dat <- na.omit(entr.dat)


# >>> GO analysis >>>
GO <- enrichGO(entr.dat$ENTREZID,  # GO富集分析
               OrgDb = GO_database,
               keyType = "ENTREZID",  # 设定读取的gene ID类型
               ont = "ALL",  # (ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
               pvalueCutoff = 1,  # 设定p值阈值
               qvalueCutoff = 0.05,  # 设定q值阈值
               readable = T)
KEGG <- enrichKEGG(entr.dat$ENTREZID,  # KEGG富集分析
                   organism = KEGG_database,
                   pvalueCutoff = 1,
                   qvalueCutoff = 0.05)

# >>> draw enrich >>>
go.bar <- barplot(GO, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free")#柱状图
kegg.bar <- barplot(KEGG,showCategory = 40,title = 'KEGG Pathway')
go.dot <- dotplot(GO, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free")#点状图
kegg.dot <- dotplot(KEGG)
# >>> save plot >>>
for (form in output.fileformat) {
  ggsave(filename = paste0(output.dir, "go_barplot", form), 
         plot = go.bar, 
         width = 8, 
         height = 10)
  ggsave(filename = paste0(output.dir, "kegg_barplot", form), 
         plot = kegg.bar, 
         width = 8, 
         height = 10)
  ggsave(filename = paste0(output.dir, "go_dotplot", form), 
         plot = go.dot, 
         width = 8, 
         height = 10)
  ggsave(filename = paste0(output.dir, "kegg_dotplot", form), 
         plot = kegg.dot, 
         width = 8, 
         height = 10)
}

# >>> network >>>
go.cnet <- enrichplot::cnetplot(GO,
                                circular = F,
                                color.params = list(edge = T),
                                node_label = "category")  # 基因-通路关联网络图
kegg.cnet <- enrichplot::cnetplot(KEGG,
                                  circular = F,
                                  color.params = list(edge = T),
                                  node_label = "category")  # circluar为指定是否环化，基因过多时建议设置为FALSE
# >>> draw network >>>
for (form in output.fileformat) {
  ggsave(filename = paste0(output.dir, "go_cnet", form), 
         plot = go.cnet, 
         width = 10, 
         height = 8)
  ggsave(filename = paste0(output.dir, "kegg_cnet", form), 
         plot = kegg.cnet, 
         width = 10, 
         height = 8)
}


# >>>>> for bulk (up and down) data use >>>>>
# >>> category by bp, cc, mf >>>
ontology.break <- list(1)
gene.ontology <- GO@result$ONTOLOGY
for (i in 2: length(gene.ontology)) {
  if (gene.ontology[i] != gene.ontology[i-1]) {
    ontology.break <- append(ontology.break, i-1)
    ontology.break <- append(ontology.break, i)
  }
}
ontology.break <- append(ontology.break, length(gene.ontology))
bp.break <- ontology.break[1: 2]
cc.break <- ontology.break[3: 4]
mf.break <- ontology.break[5: 6]

# >> 提取各个ontology的前十进行绘图 >>
circ.dat <- data.frame(ID=entr.dat$SYMBOL, logFC=log2(entr.dat$FoldChange))
GOplotIn_BP <- GO[1: 10, c(2,3,10,12)]  # 提取GO富集BP的前10行,提取ID,Description,p.adjust,GeneID四列
GOplotIn_CC <- GO[cc.break[[1]]: (cc.break[[1]]+10), c(2,3,10,12)]  # 提取GO富集CC的前10行,提取ID,Description,p.adjust,GeneID四列
GOplotIn_MF <- GO[mf.break[[1]]: (mf.break[[1]]+10), c(2,3,10,12)]  # 提取GO富集MF的前10行,提取ID,Description,p.adjust,GeneID四列
GOplotIn_BP$geneID <- str_replace_all(GOplotIn_BP$geneID, '/', ',') # 把GeneID列中的’/’替换成‘,’
GOplotIn_CC$geneID <- str_replace_all(GOplotIn_CC$geneID, '/', ',')
GOplotIn_MF$geneID <- str_replace_all(GOplotIn_MF$geneID, '/', ',')
names(GOplotIn_BP) <- c('ID','Term','adj_pval','Genes')  # 修改列名,后面弦图绘制的时候需要这样的格式
names(GOplotIn_CC) <- c('ID','Term','adj_pval','Genes')
names(GOplotIn_MF) <- c('ID','Term','adj_pval','Genes')
GOplotIn_BP$Category <- "BP"  # 分类信息
GOplotIn_CC$Category <- "CC"
GOplotIn_MF$Category <- "MF"
GOplotIn <- rbind(GOplotIn_BP, GOplotIn_CC, GOplotIn_MF)
circ_BP <- GOplot::circle_dat(GOplotIn_BP, circ.dat)  # GOplot导入数据格式整理
circ_CC <- GOplot::circle_dat(GOplotIn_CC, circ.dat) 
circ_MF <- GOplot::circle_dat(GOplotIn_MF, circ.dat)
circ_all <- GOplot::circle_dat(GOplotIn, circ.dat)
# >>> 弦表图 >>>
bp.circ <- GOCircle(circ_BP)
cc.circ <- GOCircle(circ_CC)
mf.circ <- GOCircle(circ_MF)
all.circ <- GOCircle(circ_all, 
                     nsub = 30, 
                     lfc.col = c("blue", "red"))

# >>> draw circ map >>>
for (form in output.fileformat) {
  ggsave(filename = paste0(output.dir, "bp_circ", form), 
         plot = bp.circ, 
         width = 12, 
         height = 8)
  ggsave(filename = paste0(output.dir, "cc_circ", form), 
         plot = cc.circ, 
         width = 12, 
         height = 8)
  ggsave(filename = paste0(output.dir, "mf_circ", form), 
         plot = mf.circ, 
         width = 12, 
         height = 8)
  ggsave(filename = paste0(output.dir, "all_circ", form), 
         plot = all.circ, 
         width = 40, 
         height = 30)
}

