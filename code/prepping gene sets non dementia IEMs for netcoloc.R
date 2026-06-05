library(tidyverse)
library(readxl)

# Tony prepped a sheet of genes considered to cause IEMs without neurodegen.
# here i will cleanup these gene sets in prep for running netcoloc.

# obtain the genes which are considered to cause childhood dementia.
CD_all = list.files(
  "data/genesets/v2/adultchild",
  pattern = ".csv",
  recursive = T,
  full.names = TRUE
) %>%
  sapply(function(x) {
    x %>%
      read_delim(delim = ",") %>%
      pull(1)
  }, simplify = F) %>%
  setNames(names(.) %>%
             str_remove(pattern = "data/genesets/.+/") %>%
             str_remove(pattern = ".csv$") %>%
             str_remove(pattern = "_v2")
  ) %>%
  .$CD_all


# import the gene lists from tony. tidy up ready for netcoloc.
nonCD_IEMs <- excel_sheets("data/genesets/v2/ClinGen_IEM.xlsx") %>%
  sapply(function(x)
  {
    read_xlsx(
      "data/genesets/v2/ClinGen_IEM.xlsx",
      sheet = x
    )
  },
  simplify = F) %>%
  lapply(function(x)
  {
  x %>%
      dplyr::filter(
        Classification %in% c("Definitive", "Strong"), # only retain strong evidence
        !(.[[2]] %in% CD_all) # omit genes which cause CD according to elvigdeg et al.
        ) %>%
      dplyr::select(GENE = 2) %>%
      unique %>%
      na.omit
  })

nonCD_IEMs$`FA Oxidation Disorders` <- read_xlsx(
  "data/genesets/v2/ClinGen_IEM.xlsx",
  sheet = "FA Oxidation Disorders"
) %>%
  dplyr::filter(
    Classification %in% c("Definitive", "Strong"), # only retain strong evidence
    !(Gene %in% CD_all) # omit genes which cause CD according to elvigdeg et al.
  ) %>%
  dplyr::select(GENE = 1) %>%
  unique %>%
  na.omit




# write out as an R object
saveRDS(
  nonCD_IEMs,
  "data/R_objects/genesets_nonCD_IEMs.rds"
)

# write out as csv files - individual
names(nonCD_IEMs) %>%
  lapply(function(x) {
    write.csv(
      nonCD_IEMs[[x]],
      paste0("data/genesets/v2/",x, ".csv"),
      row.names = F
    )
      })

# concat version
nonCD_IEMs %>%
  unlist %>%
  unname %>%
  as_tibble() %>%
  dplyr::select(GENE = 1) %>%
  write.csv(
    "data/genesets/v2/non_CD_IEMs_all.csv",
    row.names = F
  )




