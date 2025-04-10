#50 um
library(dplyr)

output.dir <- "./res/Norm/"

if (!dir.exists(output.dir)) {
  dir.create(output.dir, recursive = TRUE)  # recursive = TRUE 会创建所有父目录
  cat("目录已创建：", output.dir, "\n")
} else {
  cat("目录已存在：", output.dir, "\n")
}

gene.dat.path <- "./data/origin_data.txt"
sample.info.path <- "./data/sample_info.txt"

################# 改颜色，几组就改几个 ###############################
# 格式为rep(颜色, 每一组的重复实验次数)
# mycol <- c(rep("#ed6437", 3), rep("#f7b05a", 3), 
#            rep("#a6cc36", 3), rep("#b0d287", 3)) 

# mycol <- c(rep("#ed6437", 3), rep("#f7b05a", 3), 
#            rep("#a6cc36", 3))

mycol <- c(rep("#ed6437", 5), rep("#a6cc36", 5))



filter.threshold <- 0.6  # threshold that determine whether to imputate

# >> Normalize and Imputation >>
################## 里面的参数可能要改 ################################
source("./Normalization.R", local = T)



gene.dat.na <- read.delim(gene.dat.path, check.names = F)
sample.info <- read.delim(sample.info.path)
Normalization(output.dir, 
              gene.dat.na, 
              sample.info, 
              mycol, 
              filter.threshold)


# >>> PCA config >>>
output.dir <- "./res/PCA/"
gene.imputed.path <- "./res/Norm/gene_data_imputation.txt"

# >> PCA >>
gene.dat.imputed <- read.delim(gene.imputed.path, check.names = F)
source("./PCA.R", local = T)
PCA(gene.dat.imputed, sample.info, output.dir)


source("./p.R", local = T)
output.dir <- "./res/"
###################### 改 ##########################################
###################### 分成多个dataframe 分别做 ####################
control.group.name <- "CTRL"  # control叫什么就改成什么

# 选择 data 10um
gene.remove.batch <- gene.dat.imputed[, 1:7]
sample.info.batch <- sample.info[1:6, ]

# 选择 data 50um
gene.remove.batch <- gene.dat.imputed[, c(1:4, 8:10)]
sample.info.batch <- sample.info[c(1:3, 7:9), ]

gene.dat.p <- tTest(gene.remove.batch, sample.info.batch, control.group.name)
write.csv(gene.dat.p, 
          paste0(output.dir, "gene_data_p.csv"), 
          row.names = F)
write.table(gene.dat.p, 
            paste0(output.dir, "gene_data_p.txt"), 
            sep = '\t', 
            quote = F, 
            row.names = F)

source("./coveragePlot.R", local = T)
gene.bin <- calculateCoverage(gene.dat.p)
threshold <- thresholdInference(gene.bin)
coverage <- coveragePlot(gene.bin)
coverage
ggsave(paste0(output.dir, "Coverage.pdf"), 
       coverage, 
       width = 12, 
       height = 6)

source("./volcano.R", local = T)
gene.noms <- read.csv("./res/Norm/data_noms.csv")
colnames(gene.noms)[2] <- "GeneName"
gene.dat.sign <- regulationGrouping(gene.dat.p, 
                                    foldchange.threshold = threshold, 
                                    p = "p", 
                                    p.threshold = 0.05)
volcano <- volcanoPlot(gene.dat.sign, 
                       gene.noms = gene.noms, 
                       foldchange.threshold = threshold, 
                       p = "p", 
                       p.threshold = 0.05, 
                       show.annotation = c(10, 10))

gene.dat.sign <- merge(gene.noms[-4], gene.dat.sign, 
                       by.x = "Accession", by.y = "Gene")
write.csv(gene.dat.sign, 
          paste0(output.dir, "gene_data_sign.csv"), 
          row.names = F)
ggsave(paste0(output.dir, "volcano.pdf"), 
       volcano, 
       width = 12, 
       height = 12)


# >> GOKEGG >>
for (select_top in c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120)) {
  source("GO.R", local = T)
  source("KEGG.R", local = T)
  output.dir <- paste0("./res/GOKEGG_", select_top, "/")
  if (!dir.exists(output.dir)) {
    dir.create(output.dir, recursive = TRUE)  # recursive = TRUE 会创建所有父目录
    cat("目录已创建：", output.dir, "\n")
  } else {
    cat("目录已存在：", output.dir, "\n")
  }
  
  go.background <- read.csv("./data/all_uniprot_go_background.csv", header = T)
  kegg.background <- read.delim("./data/path.txt")
  colnames(gene.dat.sign)[1] <- "Gene"
  gene.up <- gene.dat.sign[gene.dat.sign$Group == "up", ]
  gene.up <- gene.up %>%
    arrange(desc(FC)) %>%
    slice(1:select_top)
  
  gene.down <- gene.dat.sign[gene.dat.sign$Group == "down", ]
  gene.down <- gene.down %>%
    arrange((FC)) %>%
    slice(1:select_top)
  
  gene.all <- rbind(gene.up, gene.down)
  gene.list <- list(gene.up, gene.down, gene.all)
  gene.sig <- c("up", "down", "all")
  for (i in 1: 3) {
    sig <- gene.sig[i]
    gene.dat.tmp <- gene.list[[i]]
    go.res <- runGOAnalysis(gene.dat.tmp, 
                            go.background)
    write.csv(go.res, paste0(output.dir, "go_result_", sig, ".csv"), row.names = F)
    go.barplot <- goBarPlot(go.res)
    go.barplot
    ggsave(paste0(output.dir, "GObarplot_", sig, ".pdf"), 
           go.barplot, 
           width = 12, 
           height = 8)
    ggsave(paste0(output.dir, "GObarplot_", sig, ".tiff"), 
           go.barplot, 
           width = 12, 
           height = 8)
    
    kegg.res <- runKEGGAnalysis(gene.dat.tmp, kegg.background)
    write.csv(kegg.res, paste0(output.dir, "kegg_result_", sig, ".csv"), row.names = F)
    kegg.plot <- drawKEGGEnrichment(kegg.res)
    ggsave(paste0(output.dir, "KEGGdotplot_", sig, ".pdf"), 
           kegg.plot$dot, 
           width = 12, 
           height = 12)
    ggsave(paste0(output.dir, "KEGGdotplot_", sig, ".tiff"), 
           kegg.plot$dot, 
           width = 12, 
           height = 12)
  }
}

library(dplyr)
gsea.output.dir <- "./res/GSEA/"

if (!dir.exists(gsea.output.dir)) {
  dir.create(gsea.output.dir, recursive = TRUE)  # recursive = TRUE 会创建所有父目录
  cat("目录已创建：", gsea.output.dir, "\n")
} else {
  cat("目录已存在：", gsea.output.dir, "\n")
}

term2gene <- kegg.background %>% select(PATH, UNIPROT)
term2name <- kegg.background %>% select(PATH, PATHNAME)
source("./GSEA.R", local = T)
gsea.go <- runGSEA(gene.dat.p, 
                   term2gene = term2gene, 
                   term2name = term2name)
write.csv(gsea.go, paste0(gsea.output.dir, "gsea_result_go.csv"), row.names = F)

gseaPlot(gsea.go, gsea.output.dir)


