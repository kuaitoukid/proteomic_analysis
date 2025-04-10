library(stats)
library(reshape2)
library(dplyr)
library(stringr)

# >>> read data >>>
gene.dat <- read.delim("./data/1_gene_data.txt")
gene.dat <- read.delim("./data/origin_data.txt")
gene.dat <- read.delim("../DIA analysis code/data/origin_data.txt")

# >>> quality control >>>
# gene.dat <- gene.dat[which(gene.dat$Sum.PEP.Score >= 5), ]
# gene.dat <- subset(gene.dat, select = -Sum.PEP.Score)
gene.dat <- na.omit(gene.dat)

# >>> calculate p value and FC >>>
stats.dat <- gene.dat[, -1] %>%
  rowwise() %>% 
  mutate(
    FC_BA = mean(
      c_across(starts_with("E")))/mean(c_across(starts_with("C"))), 
    p_AB = t.test(
      c_across(starts_with("C")), c_across(starts_with("E")))$p.value
  ) %>%
  mutate(
    FC_CA = mean(
      c_across(starts_with("E")))/mean(c_across(starts_with("C"))), 
    p_AC = t.test(
    c_across(starts_with("C")), c_across(starts_with("E")))$p.value
  )

stats.dat$p.adjust_AB <- p.adjust(stats.dat$p_AB, method = "bonferroni")
# stats.dat$p.adjust_AC <- p.adjust(stats.dat$p_AC, method = "bonferroni")
stats.dat$q_AB <- p.adjust(stats.dat$p_AB, method = "fdr")
# stats.dat$q_AC <- p.adjust(stats.dat$p_AC, method = "fdr")

gene.dat <- data.frame(Gene = gene.dat$Accession, stats.dat)
write.csv(gene.dat, "./res/1_gene_data_p.txt")
