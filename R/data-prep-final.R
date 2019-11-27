library(tidyverse)
library(readxl)

############################################################
## Read in data and calculate some preliminary statistics ##
############################################################

# Evaluate in environment to avoid global cluttering
data_env <- new.env(parent = globalenv())
evalq({

bulk_tpm <- read_tsv(
  "data/tpmMatrix.txt", skip = 1,
  col_names = c("gene", 1:74),
  col_types = cols(
  .default = col_double(),
  gene = col_character()
  ))

bulk_tmm <- read_tsv(
  "data/tmmMatrix.txt", skip = 1,
  col_names = c("gene", 1:74),
  col_types = cols(
  .default = col_double(),
  gene = col_character()
  ))

bulk_counts <- read_tsv(
  "data/countsMatrix.txt", skip = 1,
  col_names = c("gene", 1:74),
  col_types = cols(
  .default = col_double(),
  gene = col_character()
  ))

total_counts <- bulk_counts %>%
  select(-gene) %>%
  summarise_all(sum) %>%
  as_vector() %>%
  unname()

# In the following: sc = single cell
sc_tpm <- read_tsv(
  "data/scProfiles.txt",
  col_types = cols(
  .default = col_double(),
  X1 = col_character()
  )) %>%
  dplyr::rename(gene = X1)

# Merge bulk and sc data
bulk_sc_tpm <- inner_join(bulk_tpm, sc_tpm)
bulk_sc_tpm_mat <- t(as.matrix(bulk_sc_tpm %>% select(-gene)))

# Read in the design matrix
dm_in <- read_xlsx(
  "data/DesignMatrix.xlsx",
  sheet = "DesignMatrix", range = "B3:DC15",
  col_names = c("variable", 1:105),
  col_types = c("text", rep.int("numeric", 105))) %>%
  mutate_at(vars(-variable), as.integer)

dm_mat <- t(as.matrix(dm_in[,-1]))
colnames(dm_mat) <- dm_in$variable
dm <- as_tibble(dm_mat)

is_bulk <- ifelse(dm$`Bulk=1`, "B", "SC") %>%
  factor(levels = c("B", "SC"))

cell_types <- ifelse(dm$`Cell type B=1` == 1L, "B", "T") %>%
  factor(levels = c("B", "T"))

labs <- dm$`Lab`

sub_cell_types <- dm$`Sub Cell type`

tissues <- dm$`Tissue`

sample_names <- colnames(bulk_sc_tpm %>% select(-gene))

}, env = data_env)
