# Checking MPS III zebrafish netcoloc.
library(homologene)
library()

toptab_ac_cqn <- readRDS("~/analyses/2022_MPSIII_3mBrainRNAseq-main/data/R_objects/toptab_ac_cqn.rds")
B2_toptab_cqn <- readRDS("~/analyses/2022_MPSIII_3mBrainRNAseq-main/data/R_objects/B2_toptab_cqn.rds")

DE.zf = B2_toptab_cqn$`MPS-IIIB` %>%
  dplyr::filter(DE == T) %>%
  .$gene_name

DE.zf.converted2human = homologene(DE.zf, inTax = "7955", outTax = "9606") %>%
  .$`9606`

lysoDLB$gene_name %in% DE.zf.converted2human %>% summary

lysoPD$gene_name %in% DE.zf.converted2human %>% summary

logcpm <- readRDS("~/analyses/2022_MPSIII_3mBrainRNAseq-main/data/R_objects/logcpm.rds")


