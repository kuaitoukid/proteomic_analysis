library(dplyr)
library(reshape2)
library(car)
library(tidyverse)
library(stringr)
library(tidyr)

# >>> calculate p and FC >>>
#' @param gene.dat gene data, filtered by qualityControl and separateNominal
#' @details
#' gene dat should and should only have columns: Accession, ...
#' @param sample.info sample info, Sample colnames must as same as abundance rownames
#' @param control.group.name Group name of Control
tTest <- function(gene.dat, sample.info, control.group.name) {
  gene.dat.melt <- melt(gene.dat, 
                        id.vars = "Accession", 
                        variable.name = "Group", 
                        value.name = "Abundance")
  gene.dat.melt <- merge(gene.dat.melt, sample.info, 
                         by.x = "Group", by.y = "Sample", 
                         suffixes = c("", "new"))
  groups <- unique(sample.info$Group)
  
  gene.dat.p.melt <- gene.dat.melt %>% 
    group_by(Accession) %>%
    mutate(Levene.p = suppressWarnings(leveneTest(Abundance ~ Groupnew, center = mean))$`Pr(>F)`[1]) %>% 
    mutate(p = if_else(Levene.p >= 0.05, 
                       t.test(Abundance ~ Groupnew, var.equal = T)$p.value, 
                       t.test(Abundance ~ Groupnew, var.equal = F)$p.value)) %>%
    mutate(FC = mean(Abundance[Groupnew == setdiff(groups, control.group.name)]) / mean(Abundance[Groupnew == control.group.name])) %>% 
    ungroup()
  
  gene.dat.p <- gene.dat.p.melt %>% 
    pivot_longer(-c(Accession, Group, Groupnew)) %>% 
    mutate(name.new = if_else(.$name == "Abundance", paste0("Abundance_", Group), .$name)) %>% 
    select(-c(Groupnew, name, Group)) %>% 
    pivot_wider(names_from = name.new, 
                values_from = value, 
                values_fn = mean) %>% 
    rename_with(~ str_remove(., "Abundance_"), starts_with("Abundance_")) %>% 
    select(c("Accession", sample.info$Sample, "FC", "p", "Levene.p"))
  
  # >>> multiple testing >>>
  gene.dat.p$p.adjust <- p.adjust(gene.dat.p$p, method = "bonferroni") 
  gene.dat.p$q <- p.adjust(gene.dat.p$p, method = "fdr")
  colnames(gene.dat.p)[1] <- "Gene"
  return(gene.dat.p)
}