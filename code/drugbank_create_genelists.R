# messing around with the drugbank annotation

# library(drugbankR)
#
# # read in the drugbank db and convver to ta df
# # This process may take about 20 minutes.
# drugbank_dataframe <- dbxml2df(xmlfile="random/drugbank/full database.xml", version="5.1.16")
#
# drugbank_dataframe %>% colnames
#
# drugbank_dataframe$groups %>% head(20)
#
# This is isnt very easy to use.... try another way

#install.packages("dbparser")
library(dbparser)
library(tidyverse)
library(magrittr)

library(qs)

dbparser::drug_node_options()
dbparser::cett_nodes_options()

full = parseDrugBank(
  "random/drugbank/full database.xml",
  # specifiy what i want
  drug_options = dbparser::drug_node_options(),
  cett_options = dbparser::cett_nodes_options()
)

full$drugs$pathway$pathway_drugs
full$drugs$pathway$general_information$name %>% unique %>%
  grep(pattern = "auto", ignore.case = T, value = T)

full$drugs$pathway$general_information %>%
  dplyr::filter(
    grepl(name, pattern = "auto", ignore.case = T)
  )

full$cett$targets$polypeptides$go %>%
  dplyr::filter(
    grepl(description, pattern = "lyso")
    )


# plot the number of drugs per ATC code
full$drugs$atc_codes %>%
  ggplot(
    aes(x = level_4)
  ) +
  geom_bar() +
  coord_flip()

# obtain drugs annotated to be in the lysoosme by:
# pathway name
lyso.drugs.pathway.general.information <-
  full$drugs$pathway$general_information %>%
  dplyr::filter(
    grepl(name, pattern = "lysosom|autophagy|mTOR|TFEB|sphingolipid|glycosaminoglycan|mucopolysaccharid|endosom|vacuol", ignore.case = T)
  ) %>%
  .$drugbank_id %>%
  unique

## MoA
lyso.drugs.pharmacology.MoA <- full$drugs$pharmacology %>%
  dplyr::filter(
    grepl(mechanism_of_action, pattern = "lysosom|autophagy|mTOR|TFEB|sphingolipid|glycosaminoglycan|mucopolysaccharid|endosom|vacuol", ignore.case = T)
  ) %>%
  .$drugbank_id %>%
  unique

## GO of target
bp_terms <- c(
  "lysosom", "autophagy", "autophagic", "vacuo.+acidification",
  "endolysosomal", "phagolysosome", "lysosomal exocytosis",
  "lysosomal biogenesis", "lysosomal storage"
)

lyso.targets = full$cett$targets$polypeptides$go %>%
  dplyr::filter(
    grepl(description, pattern = paste(bp_terms, collapse = "|"),
                                        ignore.case = TRUE)
    ) %>%
  pull(target_id)

lyso.drugs.targets = full$cett$targets$general_information %>%
  dplyr::filter(target_id %in% lyso.targets) %>%
  pull(drugbank_id) %>%
  unique





#all_lyso_drugs <-
lyso.drugs = c(lyso.drugs.pathway.general.information,   # GO cellular component, targets
    lyso.drugs.pharmacology.MoA,  # subcellular location text, targets
    lyso.drugs.targets     # routes A+B across all CETT nodes
    # drug_ids_route_d,   # GO biological process
    # lyso_pathways,      # pathway membership
    # lyso_freetext,      # free-text mechanism/description fields
    # drug_ids_route_g,   # ATC codes
    # lyso_disease_drugs  # disease/category annotations
  ) %>% unique() %>% length



full %>%
  saveRDS("random/drugbank/fulldb.rds")
# drugbank = drugbank_db.atc.class.categ.targets$cett$targets$general_information  %>%
#   dplyr::select(target_id, name, drugbank_id) %>%
#   left_join(drugbank_db.atc.class.categ.targets$cett$targets$polypeptides$external_identy) %>%
#   dplyr::filter(resource == "GeneCards") %>%
#   left_join(drugbank_db.atc.class.categ.targets$drugs$general_information, by = "drugbank_id") %>%
#   left_join(drugbank_db.atc.class.categ.targets$drugs$atc_codes, by = "drugbank_id")
#
#
# drugbank %>%
#   split(
#     f = .$level_4
#   ) %>%
#   lapply(pull, "identifier")
#
#
#   drugbank %>%
#     saveRDS("random/drugbank/joined_atc_df.rds")


