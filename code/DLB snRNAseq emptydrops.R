# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Processing empty drops
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# I want to run empty drops through a minimal renv. so this is done here.

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# pkgs
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

library(tidyverse)
library(magrittr)
library(GEOquery)
library(DropletUtils)
library(Matrix)
library(BiocParallel)
library(qs)

# read in data ------------------------------------------------------------

# read in smaple metadata
gse = getGEO(filename="data/public data/GSE178146_RAW/GSE178146_series_matrix.txt")

meta <- pData(gse) %>%
  as_tibble %>%
  dplyr::select(
    sample = title,
    geo_accession,
    region = source_name_ch1,
    disease = characteristics_ch1
  ) %>%
  mutate(disease = str_remove(disease, pattern = "disease state: "))

# get the paths into a named vector ready for reading into emptydrops.
cellranger.parent.dir = "data/public data/GSE178146_RAW/cellranger"

sample.paths.raw =
  grep(list.dirs(cellranger.parent.dir), pattern = "raw", value = TRUE)

sample.ids <- sample.paths.raw %>%
  str_remove(pattern = "_raw_feature_bc_matrix") %>%
  str_remove(pattern = "data/public data/GSE178146_RAW/cellranger/")

sample.paths.raw <- setNames(sample.paths.raw, sample.ids)

# read in raw counts
# my laptop will run this using lapply
raw.data.counts <- lapply(
  sample.paths.raw,
  function(matrixPath) {
    cat("Loading from:", matrixPath, "\n")

    # Read count matrix
    x <- as(readMM(file.path(matrixPath, 'matrix.mtx.gz')), 'dgCMatrix')

    # Feature names
    genes <- read.delim(file.path(matrixPath, 'features.tsv.gz'), header = FALSE)
    rownames(x) <- genes[, 2]

    # Cell barcodes
    barcodes <- read.delim(file.path(matrixPath, 'barcodes.tsv.gz'), header = FALSE)
    colnames(x) <- barcodes[, 1]

    # Prefix colnames with sample ID (basename or extracted from path)
    sample_id <- basename(matrixPath)
    colnames(x) <- paste(sample_id, colnames(x), sep = "_")

    return(x)
  }
  )

# emptydrops --------------------------------------------------------------

# emptydrops
# original way i ran it on 14 samples
# this way isnt the best, since it opens 4 complete R sessions - too much memory
# set.seed(1234)
# ed.out <- raw.data.counts %>%
#   mclapply(
#     function(x) {
#       x %>%
#         emptyDrops(lower=300, retain=300, )
#     },
#     mc.cores = 4
#   )

# Will use bioc parallele instead.
# a parallel parameter object using SnowParam for stability
param <- SnowParam(workers = 4, type = "SOCK")

# Use bplapply (BiocParallel's parallel apply)
# this only took like 10 mins! amazing
ed.out <- bplapply(
  raw.data.counts,
  function(x) {
    # set a seed for reproducibility
    set.seed(1234)
    # run emptydrops
    x %>%
      emptyDrops(lower = 300, retain = 300)
  },
  BPPARAM = param
)

E.keep <- ed.out %>%
  lapply(function(x)
  {x$FDR <= 0.001}
  )

E.keep.na = E.keep

E.keep.na <- lapply(E.keep.na, function(x) {
  x[is.na(x)] <- FALSE
  return(x)
})

summary_lengths <- lapply(E.keep.na, summary)

summary_lengths %>%
  bind_rows(.id = "sample") %>%
  dplyr::select(
    sample, numCells = `TRUE`
  ) %>%
  mutate(numCells = as.numeric(numCells)) %>%
  left_join(meta) %>%
  ggplot(
    aes(x = sample, y = numCells, fill = disease)
  ) +
  geom_col() +
  facet_wrap(~disease, scales = "free_x", nrow = 1) +
  labs(
    y = "# of cells",
    title = "Number of droplets after filtering with emptyDrops()"
  )

# perform the filtering
filtered.matrices <- map2(raw.data.counts, E.keep.na, function(mat, keep) {
  keep[is.na(keep)] <- FALSE
  mat[, keep]
})


# Create Seurat object
library(Seurat)
data.seurat <- CreateSeuratObject(
  counts = filtered.matrices,
  min.cells = 5,
  project = "DLB_PD_filtered"
)

colnames(data.seurat) %<>%
  str_remove("_raw_feature_bc_matrix")

# save RDS for importing into rmarkwon later
data.seurat %>%
  qsave("data/public data/GSE178146_RAW/R_objects/data_seurat_DLB_PD_snRNAseq_emptydrops.qs",
      preset = "fast",
      nthreads = 4)
