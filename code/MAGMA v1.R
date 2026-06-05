
files = list.files("data/MAGMA_v2", full.names = T)


# MAGMA top model ---------------------------------------------------------
# Input file:
# UPDATED_2026_CD_genes_MAGMA_TOP_conservative_boundaries_across_all_pheno.txt

## Our primary gene-based test, forms basis for supplementary figure with boxplot and heatmap
CD_merged = files[grep(files,
           pattern = "MAGMA_TOP_conservative_boundaries_across_all_pheno.txt")] %>%
  read.delim() %>%
  as_tibble()

## Count how many genes surpass either threshold
## LBD = 1
CD_merged %>% filter(DLB < 5.128205e-05)
## LOAD = 1
CD_merged %>% filter(AD < 5.128205e-05)
## FTD = none
CD_merged %>% filter(FTLD < 5.128205e-05)
## ALS = 2
CD_merged %>% filter(ALS < 5.128205e-05)
## PD = 8
CD_merged %>% filter(PD < 5.128205e-05)

Boxplot <- CD_merged %>%
  gather(key = "disease", value = "p", -starts_with("GENE"))

p1 = ggplot(Boxplot, aes(x=disease, y=-log10(p), fill=disease)) +
  geom_boxplot(alpha = 0.7) +
  geom_hline(yintercept = 4.290035, lty = "dashed") +
  geom_hline(yintercept = 6.250152, lty = "dashed") +
  theme_bw() +
  scale_fill_viridis(discrete = T) +
  xlab(" ") + ylab("-log10 MAGMA (top) P") +
 # ggtitle("Gene based association of childhood dementia genes") +
  theme(legend.position = "none",
         text = element_text(size = 12))

## Make a heat-map of -log10(P) of all Bonferroni sig genes
Bonferroni <- CD_merged %>% filter(PD < 5.128205e-05 | ALS < 5.128205e-05 | DLB < 5.128205e-05 | AD < 5.128205e-05)

Merged_mat <- as.matrix(Bonferroni)

rownames(Merged_mat) <- Merged_mat[, 1]

Merged_mat <- Merged_mat[, -1]


Merged_mat <- `dimnames<-`(`dim<-`(as.numeric(Merged_mat), dim(Merged_mat)), dimnames(Merged_mat))
Merged_mat <- -log10(Merged_mat)

col_fun = colorRamp2(c(0, 10), c("white", "purple"))


HT <- Heatmap(
  Merged_mat,
  name = "neglog10P",
  rect_gp = gpar(col = "white", lwd = 2),
  show_column_dend = TRUE, show_row_dend = TRUE,
  clustering_distance_rows = "pearson",

  row_dend_width = unit(1, "cm"),
  col = col_fun,
  row_names_gp = grid::gpar(fontsize = 12),
  column_names_gp = grid::gpar(fontsize = 12, hjust = 0),
  column_names_centered = T,
  column_names_rot = 360,
  column_title_gp = grid::gpar(fontsize = 11.5, fontface="bold"),
  cell_fun = function(j, i, x, y, w, h, fill) {
    if(Merged_mat[i, j] > 4.290035 && Merged_mat[i, j] < 6.250152) {
      grid.text("*", x, y)
    }
    if(Merged_mat[i, j] > 6.250152) {
      grid.text("**", x, y)
    }
  }
)

ht_grob <- grid.grabExpr({
  ComplexHeatmap::draw(HT)
})
png("output/plots4paperv2/suppMagma_top.png", width = 20, height = 10, units = "cm", res = 300)
cowplot::plot_grid(p1, ht_grob, ncol = 2,
                   labels = c("A", "B"), label_size = 20, rel_widths = c(0.4, 0.6))
dev.off()

# MAGMA multi model  ------------------------------------------------
# Our alternative gene-based test, forms basis for second supplementary figure (MAGMA multi) with boxplot and heatmap

# UPDATED_2026_CD_genes_MAGMA_multi_conservative_boundaries_across_all_pheno.txt

CD_merged.multi = files[grep(files,
                       pattern = "MAGMA_multi_conservative_boundaries_across_all_pheno.txt")] %>%
  read.delim() %>%
  as_tibble()

## Count how many genes surpass either threshold
## LBD = 0
CD_merged.multi %>% filter(DLB < 5.128205e-05)
## LOAD = 0
CD_merged.multi %>% filter(AD < 5.128205e-05)
## FTD = none
CD_merged.multi %>% filter(FTLD < 5.128205e-05)
## ALS = 3
CD_merged.multi %>% filter(ALS < 5.128205e-05)
## PD = 8
CD_merged.multi %>% filter(PD < 5.128205e-05)

Boxplot.multi <- CD_merged.multi %>%
  gather(key = "disease", value = "p", -starts_with("GENE"))

p2 = ggplot(Boxplot.multi, aes(x=disease, y=-log10(p), fill=disease)) +
  geom_boxplot(alpha = 0.7) +
  geom_hline(yintercept = 4.290035, lty = "dashed") +
  geom_hline(yintercept = 6.250152, lty = "dashed") +
  theme_bw() +
  scale_fill_viridis(discrete = T) +
  xlab(" ") + ylab("-log10 MAGMA (top) P") +
  # ggtitle("Gene based association of childhood dementia genes") +
  theme(legend.position = "none",
        text = element_text(size = 12))

## Make a heat-map of -log10(P) of all Bonferroni sig genes
Bonferroni.multi <- CD_merged.multi %>%
  filter(PD < 5.128205e-05 | ALS < 5.128205e-05 | DLB < 5.128205e-05 | AD < 5.128205e-05)

Merged_mat.multi <- as.matrix(Bonferroni.multi)

rownames(Merged_mat.multi) <- Merged_mat.multi[, 1]

Merged_mat.multi <- Merged_mat.multi[, -1]


Merged_mat.multi <- `dimnames<-`(`dim<-`(as.numeric(Merged_mat.multi), dim(Merged_mat.multi)), dimnames(Merged_mat.multi))
Merged_mat.multi <- -log10(Merged_mat.multi)

col_fun = colorRamp2(c(0, 10), c("white", "purple"))


HT.multi <- Heatmap(
  Merged_mat.multi,
  name = "neglog10P",
  rect_gp = gpar(col = "white", lwd = 2),
  show_column_dend = TRUE, show_row_dend = TRUE,
  clustering_distance_rows = "pearson",

  row_dend_width = unit(1, "cm"),
  col = col_fun,
  row_names_gp = grid::gpar(fontsize = 12),
  column_names_gp = grid::gpar(fontsize = 12, hjust = 0),
  column_names_centered = T,
  column_names_rot = 360,
  column_title_gp = grid::gpar(fontsize = 11.5, fontface="bold"),
  cell_fun = function(j, i, x, y, w, h, fill) {
    if(Merged_mat.multi[i, j] > 4.290035 && Merged_mat.multi[i, j] < 6.250152) {
      grid.text("*", x, y)
    }
    if(Merged_mat.multi[i, j] > 6.250152) {
      grid.text("**", x, y)
    }
  }
)

ht_grob2 <- grid.grabExpr({
  ComplexHeatmap::draw(HT.multi)
})

png("output/plots4paperv2/suppMagma_multi.png", width = 20, height = 10, units = "cm", res = 300)
cowplot::plot_grid(p2, ht_grob2, ncol = 2,
                   labels = c("A", "B"), label_size = 20, rel_widths = c(0.4, 0.6))
dev.off()



# MAGMA GSA ---------------------------------------------------------------

# top SNP model -----------------------------------------------------------

# input file:
# UPDATED_2026_MAGMA_top_GSA_FDR_corr.txt

# Input for supplementary figure with points for -log10 P per disorder (primary MAGMA method)

p.topSNP.GSA = files %>%
  grep(pattern = "top_GSA", value = T) %>%
  read.delim() %>%
  as_tibble %>%
  mutate(
    type = case_when(
      grepl(Set, pattern = "propagated") ~ "Propagated network",
      TRUE ~ "Seed genes only"
    ),

    Disorder = case_when(
      Disorder == "LBD" ~ "DLB",
      Disorder == "FTD" ~ "FTLD",
      TRUE ~ Disorder
    ),

    Disease = str_remove(Set, pattern = " \\(propagated\\)"),
    # convert to full anames:
    Disease = case_when(

      grepl(Disease, pattern = "NBIA") ~ "Neurodegeneration with brain iron accumulation",

      grepl(Disease, pattern = "All genes") ~ "All childhood dementia genes",

      grepl(Disease, pattern = "Amino ac") ~ "Disorders of amino acid & other organic acid metabolism",

      grepl(Disease, pattern = "Leukodystrophies") ~ "Leukodystrophies not otherwise categorised",

      grepl(Disease, pattern = "Lysosomal") ~ "Lysosomal disorders",

      grepl(Disease, pattern = "Mitochondrial") ~ "Mitochondrial disorders",

      grepl(Disease, pattern = "Vitamin") ~ "Vitamin-responsive inborn errors of metabolism",

      grepl(Disease, pattern = "Peroxisomal") ~ "Peroxisomal disorders",

      grepl(Disease, pattern = "ther lipid") ~ "Other disorders of lipid metabolism & transport",

      grepl(Disease, pattern = "Other IEM") ~ "Other inborn errors of metabolism",

      grepl(Disease, pattern = "Mineral absorption") ~ "Disorders of mineral absorption & transport",

      grepl(Disease, pattern = "Not categorised") ~ "Diseases not otherwise categorised"

    ) %>% str_wrap(width = 30)
    ) %>%
  arrange(Disease) %>%
  mutate(
    Disease = fct_inorder(Disease) %>% fct_rev(),
    type = fct_rev(type)
  ) %>%
  ggplot(
    aes(y = Disease, x = -log10(pvalue), colour = Disorder)
  ) +
  geom_point(
    size = 5,
    alpha = 0.7,
  ) +
  geom_point(
    size = 5,
    alpha = 0.7,
    shape = 1,
    colour = "grey10"
  ) +
  theme_bw() +
  facet_wrap(~type) +
  scale_colour_viridis_d(option = "viridis") +
  geom_vline(xintercept = 1.3, lty="dashed", colour = "grey50") +
  geom_vline(xintercept = 3.35, lty="dashed", colour = "red") +
  labs(
    y = NULL,
    x = "-log10(P) - MAGMA GSA (top SNP model)",
    colour = "Adult-onset NDD"
    )


# multi -------------------------------------------------------------------

# 4) UPDATED_2026_MAGMA_multi_GSA_FDR_corr.txt
#
# Input for supplementary figure with points for -log10 P per disorder (secondary MAGMA multi  method)

p.multi.GSA = files %>%
  grep(pattern = "multi_GSA", value = T) %>%
  read.delim() %>%
  as_tibble %>%
  mutate(
    type = case_when(
      grepl(Set, pattern = "propagated") ~ "Propagated network",
      TRUE ~ "Seed genes only"
    ),

    Disorder = case_when(
      Disorder == "LBD" ~ "DLB",
      Disorder == "FTD" ~ "FTLD",
      TRUE ~ Disorder
      ),

    Disease = str_remove(Set, pattern = " \\(propagated\\)"),
    # convert to full anames:
    Disease = case_when(

      grepl(Disease, pattern = "NBIA") ~ "Neurodegeneration with brain iron accumulation",

      grepl(Disease, pattern = "All genes") ~ "All childhood dementia genes",

      grepl(Disease, pattern = "Amino ac") ~ "Disorders of amino acid & other organic acid metabolism",

      grepl(Disease, pattern = "Leukodystrophies") ~ "Leukodystrophies not otherwise categorised",

      grepl(Disease, pattern = "Lysosomal") ~ "Lysosomal disorders",

      grepl(Disease, pattern = "Mitochondrial") ~ "Mitochondrial disorders",

      grepl(Disease, pattern = "Vitamin") ~ "Vitamin-responsive inborn errors of metabolism",

      grepl(Disease, pattern = "Peroxisomal") ~ "Peroxisomal disorders",

      grepl(Disease, pattern = "ther lipid") ~ "Other disorders of lipid metabolism & transport",

      grepl(Disease, pattern = "Other IEM") ~ "Other inborn errors of metabolism",

      grepl(Disease, pattern = "Mineral absorption") ~ "Disorders of mineral absorption & transport",

      grepl(Disease, pattern = "Not categorised") ~ "Diseases not otherwise categorised"

    ) %>% str_wrap(width = 30),

  ) %>%
  arrange(Disease) %>%
  mutate(
    Disease = fct_inorder(Disease) %>% fct_rev(),
    type = fct_rev(type)
  ) %>%
  ggplot(
    aes(y = Disease, x = -log10(pvalue), colour = Disorder)
  ) +
  geom_point(
    size = 5,
    alpha = 0.7,
  ) +
  geom_point(
    size = 5,
    alpha = 0.7,
    shape = 1,
    colour = "grey10"
  ) +
  theme_bw() +
  facet_wrap(~type) +
  scale_colour_viridis_d(option = "viridis") +
  geom_vline(xintercept = 1.3, lty="dashed", colour = "grey50") +
  geom_vline(xintercept = 3.35, lty="dashed", colour = "red") +
  labs(
    y = NULL,
    x = "-log10(P) - MAGMA GSA (multi SNP model)",
    colour = "Adult-onset NDD"
  )

ggarrange(p.topSNP.GSA, p.multi.GSA,
          common.legend = T, labels = "AUTO") +
ggsave(
    "output/plots4paperv2/supp/suppMagmaGSA_top&multi.png",
    width = 18, height = 10, units = "cm",
    scale = 1.5
  )
