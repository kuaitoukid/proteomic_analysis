library(clusterProfiler)
library(tidyverse)
library(enrichplot)
library(ggplot2)


#' @description
#' GSEA with self contructed database
#' @param gene.dat gene data with FC and p
#' @param term2gene PathID to UniprotID
#' @param term2name PahtID to Paht name
#' @param p should be one of p, p.adjust or q, default p
#' @param p.threshold default 0.05
runGSEA <- function(gene.dat, 
                    term2gene, 
                    term2name, 
                    p = "p", 
                    p.threshold = 0.05) {
  # >>> data filter >>>
  gene.dat <- gene.dat[gene.dat[[p]] < p.threshold, ]
  gene.dat <- subset(gene.dat, select = c(Gene, FC))
  gene.dat$FC <- log2(gene.dat$FC)
  
  # >>> sort by FC >>>
  gene.dat <- gene.dat[order(gene.dat$FC, decreasing = T), ]
  gene.list <- gene.dat$FC
  names(gene.list) <- gene.dat$Gene
  
  gsea <- GSEA(geneList = gene.list,
               TERM2GENE = term2gene,
               TERM2NAME = term2name, 
               pvalueCutoff = 1, 
               eps = 1e-50)
  
  return(gsea)
}


#' @description
#' Draw ridgeplot and gseaplot
#' @param gsea gsae result
#' @param output.dir output dir path
#' @param gsea.p.threshold default 0.05
gseaPlot <- function(gsea, 
                     output.dir, 
                     gsea.q.threshold = 0.05) {
  gsea.res <- gsea@result
  show.category <- gsea.res %>% 
    filter(qvalue < gsea.q.threshold) %>% 
    nrow()
  if (show.category > 30) {
    show.category <- 30
  }
  
  # >>> ridge plot >>>
  tryCatch({
    ridge.plot <- ridgeplot(gsea, 
                            showCategory = show.category, 
                            orderBy = "qvalue", 
                            decreasing = T) + 
      theme(text = element_text(size = 10))
    ggsave(paste0(output.dir, "ridge_plot.pdf"), 
           ridge.plot, 
           width = 8, 
           height = 8)
  }, error = function(e) {print(e)})
  
  
  gseaplot.dir <- paste0(output.dir, "gseaplot/")
  if (!dir.exists(gseaplot.dir)) {
    dir.create(gseaplot.dir)
  } else {
    warning(paste0(output.dir, "gseaplot dir exists!"))
  }
  width <- options()$width #获取显示界面行宽度
  N <- nrow(gsea.res)
  for (i in 1: nrow(gsea.res)) {
    gsea.line <- gsea.res[i, ]
    path.id <- gsea.line$ID
    path.name <- gsea.line$Description
    path.info <- paste0("ES: ", gsea.line$enrichmentScore, 
                        "\tqvalue: ", gsea.line$qvalue)
    
    gsea.plot <- gseaplot2(gsea, 
                           path.id, 
                           title = paste0(path.id, ": ", path.name, 
                                          "\n", path.info), 
                           pvalue_table = F, 
                           base_size = 6)
    ggsave(paste0(gseaplot.dir, path.id, "_", 
                  gsub(" / ", "_", path.name), ".pdf"), 
           gsea.plot, 
           width = 7, height = 4)
    ggsave(paste0(gseaplot.dir, path.id, "_", 
                  gsub(" / ", "_", path.name), ".png"), 
           gsea.plot, 
           width = 7, height = 4)
    cat('[', paste0(rep('#', i/N*width), collapse=''),
        paste0(rep('-', width - i/N*width), collapse=''),
        ']',
        round(i/N*100),'%')
    if(i==N)cat('\nDONE!\n')
    else cat('\r')
  }
}
