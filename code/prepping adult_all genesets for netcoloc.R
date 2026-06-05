## making a new adult_all.csv geneset

## KB
## 19/2/25

library(tidyverse)
library(magrittr)

# read in the genesets from will - 1 gene per GWAS hit
genesets <-
  ## ADULT ONSET
    list.files(
      path = "data/genesets/Alternate_gene_lists_for_NetColoc_senstivity_analyses",
      pattern = ".txt",
      recursive = T,
      full.names = TRUE
    ) %>%
      grep(pattern = "_list.txt", value = T) %>%
      sapply(function(x) {
        x %>%
          read_delim(delim = "t") %>%
          .$Gene
      }, simplify = F) %>%
      setNames(names(.) %>%
                 str_remove(pattern = "data/genesets/Alternate_gene_lists_for_NetColoc_senstivity_analyses/") %>%
                 str_remove(pattern = "_v2.+txt"))

genesets %>% names

tibble(
  GENE = unlist(genesets) %>%
    unname() # remove the names on the vectot
) %>%
  write_csv(
    "data/genesets/adult_all.csv"
  )




