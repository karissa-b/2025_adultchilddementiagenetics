


library(Seurat)
library(SingleCellExperiment)
library(scPower)
library(tidyverse)
library(qs)

# Data object
data.seurat.mm.Npc1 = qread("data/public data/PRJNA606815_NPC_mousebrain/R_objects/dataseurat.qs")
data.seurat.PD.DLB = qread("data/public data/GSE178146_RAW/R_objects/PD_DLB_data_seurat.qs")


# ---- Inputs ----------------------------------------------------------------
seu          <- data.seurat.PD.DLB
target_ct    <- "Neuron_Ex"   # change to whichever cell type you want
donor_col    <- "orig.ident"
min_counts   <- 3
perc_indiv   <- 0.5

# ---- Subset to target cell type --------------------------------------------
seu_ct <- subset(seu, subset = celltype == target_ct)

counts_ct <- GetAssayData(seu_ct, assay = "RNA", slot = "counts")

annot_ct <- seu_ct@meta.data %>%
  rownames_to_column("cell") %>%
  transmute(cell, cell_type = celltype, individual = .data[[donor_col]]) %>%
  as.data.frame()
rownames(annot_ct) <- annot_ct$cell

# ---- Build complete + 3 subsampled count matrices --------------------------
make_subsample <- function(prop, mat) {
  if (prop == 1) return(mat)
  target <- round(colSums(mat) * prop)
  SampleUMI(mat, max.umi = target, upsample = FALSE, verbose = FALSE)
}

count_mats <- tibble(
  matrix = c("complete", "subsampled75", "subsampled50"),
  prop   = c(1, 0.75, 0.50)
) %>%
  mutate(mat = map(prop, make_subsample, mat = counts_ct))

# ---- Step 1: count expressed genes per matrix ------------------------------
expressed_genes_df <- count_mats %>%
  mutate(
    n_cells = map_int(mat, ncol),
    pseudo  = map(mat, ~ create.pseudobulk(.x, annot_ct[colnames(.x), ],
                                           colName = "cell_type")),
    expr    = map(pseudo, calculate.gene.counts,
                  min.counts = min_counts, perc.indiv = perc_indiv),
    expressed.genes = map_int(expr, nrow)
  ) %>%
  dplyr::select(matrix, n_cells, expressed.genes)

expressed_genes_df

# ---- Step 2: NB fit per matrix ---------------------------------------------
nb_fits <- count_mats %>%
  mutate(fit = map(mat, ~ nbinom.estimation(as.matrix(.x),
                                            sizeFactorMethod = "poscounts")))

norm_mean_values <- nb_fits %>%
  mutate(means = map(fit, 1)) %>%
  select(matrix, means) %>%
  unnest(means)

disp_param <- nb_fits %>%
  mutate(disp = map(fit, 3)) %>%
  select(matrix, disp) %>%
  unnest(disp)

# ---- Step 3: gamma mixed fit per matrix ------------------------------------
gamma_fits <- count_mats %>%
  mutate(
    censored = 1 / map_int(mat, ncol),
    means    = map(matrix, ~ norm_mean_values$mean[norm_mean_values$matrix == .x]),
    gfit     = map2(means, censored,
                    ~ mixed.gamma.estimation(.x,
                                             num.genes.kept = 26000,
                                             censoredPoint  = .y))
  ) %>%
  select(matrix, gfit) %>%
  unnest(gfit)

# ---- Step 4: mean UMI per matrix -------------------------------------------
umi_values <- count_mats %>%
  transmute(matrix, mean.umi = map_dbl(mat, ~ meanUMI.calculation(as.matrix(.x))))

# ---- Step 5: parameterize gamma fits over mean UMI -------------------------
gamma_linear_fit <- gamma_fits %>%
  left_join(umi_values, by = "matrix") %>%
  convert.gamma.parameters() %>%
  umi.gamma.relation() %>%
  mutate(ct = target_ct)

# ---- Step 6: dispersion function -------------------------------------------
disp_fun <- dispersion.function.estimation(disp_param) %>%
  mutate(ct = target_ct)

# ---- Step 7: read–UMI fit --------------------------------------------------
# Approximation: scale total UMIs by an assumed UMI/read ratio (~0.4).
# Replace with real cellranger mapped reads if you have them.
approx_reads <- sum(counts_ct) / 0.4

read_umi_fit <- umi_values %>%
  left_join(count_mats %>% select(matrix, prop), by = "matrix") %>%
  mutate(transcriptome.mapped.reads = approx_reads * prop / ncol(counts_ct)) %>%
  select(matrix, mean.umi, transcriptome.mapped.reads) %>%
  umi.read.relation() %>%
  mutate(type = target_ct)

# ---- Save priors -----------------------------------------------------------
priors <- list(
  gamma.mixed.fits = gamma_linear_fit,
  disp.fun.param   = disp_fun,
  read.umi.fit     = read_umi_fit,
  expressed.genes  = expressed_genes_df
)

# ---- Example power calculation ---------------------------------------------
ct_freq <- mean(seu@meta.data$celltype == target_ct)

# power_res <- power.general.withDoublets(
#   nSamples          = 100,
#   nCells            = 1500,
#   readDepth         = 25000,
#   ct.freq           = ct_freq,
#   type              = "eqtl",
#   ref.study         = scPower::eqtl.ref.study,
#   ref.study.name    = "Blueprint (Monocytes)",
#   samplesPerLane    = 4,
#   read.umi.fit      = read_umi_fit,
#   gamma.mixed.fits  = gamma_linear_fit,
#   ct                = target_ct,
#   disp.fun.param    = disp_fun,
#   mappingEfficiency = 0.8,
#   min.UMI.counts    = min_counts,
#   perc.indiv.expr   = perc_indiv,
#   sign.threshold    = 0.05,
#   MTmethod          = "FDR"
# )
#
# print(power_res)

# ---- Sample size sweep for 80% power ---------------------------------------
# sample_sizes <- seq(2, 10, by = 2)
#
# power_sweep <- tibble(nSamples = sample_sizes) %>%
#   mutate(
#     res = map(nSamples, ~ power.general.withDoublets(
#       nSamples          = .x,
#       nCells            = 1500,
#       readDepth         = 25000,
#       ct.freq           = ct_freq,
#       type              = "de",
#       ref.study         = scPower::de.ref.study,
#       ref.study.name    =  de.ref.study$name %>% unique %>% .[2],
#       samplesPerLane    = 4,
#       read.umi.fit      = read_umi_fit,
#       gamma.mixed.fits  = gamma_linear_fit,
#       ct                = target_ct,
#       disp.fun.param    = disp_fun,
#       mappingEfficiency = 0.8,
#       min.UMI.counts    = min_counts,
#       perc.indiv.expr   = perc_indiv,
#       sign.threshold    = 0.05,
#       MTmethod          = "FDR"
#     ))
#   ) %>%
#   unnest_wider(res)

# ---- Plot ------------------------------------------------------------------
# ggplot(power_sweep, aes(nSamples, powerDetect)) +
#   geom_line() +
#   geom_point() +
#   geom_hline(yintercept = 0.90, linetype = "dashed", colour = "red") +
#   geom_vline(xintercept = n_for_80, linetype = "dashed", colour = "red") +
#   scale_y_continuous(limits = c(0, 1)) +
#   labs(x = "Number of samples",
#        y = "Detection power",
#        title = paste0("Power curve for ", target_ct,
#                       " (1500 cells, 25k read depth)")) +
#   theme_bw()


# ---- Simulate DE priors ----------------------------------------------------
# Pick effect size assumptions for your Npc1 vs WT contrast.
# foldChange is sampled from a normal on the log scale, then exponentiated.
# mean = 2 means average log2 fold change of 2 (i.e. ~4x). Adjust to taste.
n_de_genes  <- 200
fc_mean     <- 1.0   # log2 fold change mean — try 0.5, 1, 1.5 to see sensitivity
fc_sd       <- 0.5

# Gene ranks: where in the expression distribution the DE genes sit.
# Mid-to-high expression range is typical for detectable DE genes.
ranks <- uniform.ranks.interval(start = 1000, end = 7000, numGenes = n_de_genes)[1:n_de_genes]

simulated_de <- tibble(
  ranks      = ranks,
  FoldChange = effectSize.DE.simulation(mean = fc_mean,
                                        sd   = fc_sd,
                                        numGenes = n_de_genes),
  name       = "Simulated"
) %>% as.data.frame()

# ---- Sample size sweep -----------------------------------------------------
sample_sizes <- seq(2, 10, by = 2)

power_sweep <- tibble(nSamples = sample_sizes) %>%
  mutate(
    res = map(nSamples, ~ power.general.withDoublets(
      nSamples          = .x,
      nCells            = 1500,
      readDepth         = 25000,
      ct.freq           = ct_freq,
      type              = "de",
      ref.study         = simulated_de,
      ref.study.name    = "Simulated",
      samplesPerLane    = 4,
      read.umi.fit      = read_umi_fit,
      gamma.mixed.fits  = gamma_linear_fit,
      ct                = target_ct,
      disp.fun.param    = disp_fun,
      mappingEfficiency = 0.8,
      min.UMI.counts    = min_counts,
      perc.indiv.expr   = perc_indiv,
      sign.threshold    = 0.05,
      MTmethod          = "FDR"
    ))
  ) %>%
  unnest_wider(res)

ggplot(power_sweep, aes(nSamples, powerDetect)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0.80, linetype = "dashed", colour = "red") +
  geom_vline(xintercept = n_for_80, linetype = "dashed", colour = "red") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Number of samples",
       y = "Detection power",
       title = paste0("DE power for ", target_ct,
                      " (log2FC ~ N(", fc_mean, ", ", fc_sd, "))")) +
  theme_bw()
