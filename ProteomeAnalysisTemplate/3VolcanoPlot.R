library(tidyverse)
library(dplyr)

# >>> read file >>>
gene.dat <- read.csv("./res/1_gene_data_p.txt")

# >>> filter >>>
# use Group A and B as example
gene.dat$sig <- ifelse(gene.dat$FC_BA >= 1.3 & gene.dat$q_AB <= 0.05, "up", 
                  ifelse(gene.dat$FC_BA <= 0.7 & gene.dat$q_AB <= 0.05, "down", 
                         "no"
                  )
)

up.count <- nrow(gene.dat[which(gene.dat$sig == "up"), ])
down.count <- nrow(gene.dat[which(gene.dat$sig == "down"), ])
no.count <- nrow(gene.dat[which(gene.dat$sig == "no"), ])

logFC = log2(gene.dat$FC_BA)
qval = -log10(gene.dat$q_AB)
sig = gene.dat$sig
volcano.dat <- data.frame(Gene = gene.dat$Gene, 
                          logFC = logFC, 
                          qval = qval, 
                          sig = sig)

# >>> extract most far gene to sign name >>>
library(org.Hs.eg.db)
signcount.up <- 10
signcount.down <- 10
volcano.dat$distance <- volcano.dat$logFC^2 + volcano.dat$qval^2
sign.up <- volcano.dat[which(volcano.dat$sig == "up"), ]
sign.up <- head(sign.up[order(sign.up$distance, decreasing = T),], signcount.up)
sign.down <- volcano.dat[which(volcano.dat$sig == "down"), ]
sign.down <- head(sign.down[order(sign.down$distance, decreasing = T),], signcount.down)
sign <- rbind(sign.up, sign.down)
gene.map <- AnnotationDbi::select(org.Hs.eg.db, 
                                  keys=sign$Gene, 
                                  keytype="UNIPROT", 
                                  columns="SYMBOL")

colnames(gene.map)[1] <- "Gene"
sign.all <- merge(volcano.dat, 
                  gene.map,
                  by = "Gene")
sign.all <- na.omit(sign.all)

# >>> volcano plot >>>
library(ggplot2)
library(ggrepel)
library(RColorBrewer)

theme_set(theme_bw())
p <- ggplot(volcano.dat, aes(logFC, qval, color = sig)) + 
  geom_point(size = 1, aes(color = sig)) + 
  xlim(-max(abs(logFC)), max(abs(logFC))) + 
  labs(x="log2FC", y="-log10qval")
p <- p + scale_color_manual(values = c("down" = "#0072B5", "no" = 'grey', "up" = '#BC3C28'))+
  geom_hline(yintercept = c(-log10(0.05)), linetype = 4) +
  geom_vline(xintercept = c(log2(0.7), log2(1.3)), linetype = 4)
p <- p + theme(panel.grid = element_blank()) +
  theme(axis.line = element_line(size=0))
p <- p + geom_text_repel(  #  标注基因
  data = sign.all, 
  aes(label = SYMBOL),
  size = 3,
  color = "black",
  segment.color = "black", show.legend = FALSE)
p <- p + guides()
p <- p + 
  theme(axis.text = element_text(size = 10), 
  axis.title = element_text(size = 10))

p
ggsave("./res/3_volcano.pdf", p, width = 10, height = 10)