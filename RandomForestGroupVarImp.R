library(ggplot2)
library(dplyr)
library(randomForestSRC)
# The following additional packages need to be installed:
# - bootstrap
# - stringr
# - scales
# - ggpubr

data_folder <- "./data/"

# Colour-blind friendly palette
cbPalette <- c(
	"#999999", "#E69F00", "#56B4E9", "#009E73",
	"#F0E442", "#0072B2", "#D55E00", "#CC79A7")

##########################
## Setup & Convert data ##
##########################

one_hot_enc <- function(old_var, new_var_name) {
    enc <- do.call(
        cbind, lapply(sort(unique(old_var)),
                      function(i) as.integer(old_var == i)))
    colnames(enc) <- sprintf(
      "%s%d", new_var_name,
      sort(unique(old_var)))

    # Return encoded matrix as data.frame
    as.data.frame(enc)
}

# Work with covariates
# - Lab
# - Tissue
# - Bulk or SC
# - Cell Type
# - Sub Cell Type

sc_or_bulk <- readRDS(paste0(data_folder, "scOrBulk.RDS"))
celltype <- readRDS(paste0(data_folder, "cellTypes.RDS"))
subcelltype <- readRDS(paste0(data_folder, "subCellTypes.RDS"))
lab <- readRDS(paste0(data_folder, "labs.RDS"))
tissue <- readRDS(paste0(data_folder, "tissues.RDS"))

sc_or_bulk_enc <- one_hot_enc(sc_or_bulk, "sc_or_bulk")
celltype_enc <- one_hot_enc(celltype, "celltype")
subcelltype_enc <- one_hot_enc(subcelltype, "subcelltype")
lab_enc <- one_hot_enc(lab, "lab")
tissue_enc <- one_hot_enc(tissue, "tissue")

# Encoded design matrix
X <- cbind(
    sc_or_bulk_enc, celltype_enc,
    subcelltype_enc,
    lab_enc, tissue_enc)

# Load gene datasets (samples in columns)
tpm_sc_and_bulk <- readRDS(paste0(data_folder, "tpmScAndBulk.RDS"))
tmm_sc_and_bulk <- readRDS(paste0(data_folder, "tmmScAndBulk.RDS"))
quantile_sc_and_bulk <- readRDS(paste0(data_folder, "quantileScAndBulk.RDS"))

# Exclude technical replicates
# There are too few technical replicates to measure the impact of these
tech_reps <- c(3, 4, 6, 7, 11)
X <- X[-tech_reps,]
tpm_sc_and_bulk <- tpm_sc_and_bulk[,-tech_reps]
tmm_sc_and_bulk <- tmm_sc_and_bulk[,-tech_reps]
quantile_sc_and_bulk <- quantile_sc_and_bulk[,-tech_reps]
sc_or_bulk <- sc_or_bulk[-tech_reps]
celltype <- celltype[-tech_reps]
subcelltype <- subcelltype[-tech_reps]
lab <- lab[-tech_reps]
tissue <- tissue[-tech_reps]

# Plot the frequency of covariate values
tibble(
  variable = factor(colnames(X), levels = colnames(X)),
  frequency = colSums(X)) %>%
  ggplot() +
  geom_bar(aes(x = variable, y = frequency), stat = "identity") +
  labs(title = "Frequency of covariates", x = NULL, y = "Frequency") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("plots/2019-12-06-covariate-frequency.png", device = "png",
  width = 15, height = 5)

# Work in full principal component space, i.e. no information loss regarding
# the variance between samples.
# This allows us to work in 100 dimensional space instead of ~ 15000
Y <- list(
  tpm = prcomp(log2(t(tpm_sc_and_bulk) + 1))$x,
  tmm = prcomp(log2(t(tmm_sc_and_bulk) + 1))$x,
  quantile = prcomp(log2(t(quantile_sc_and_bulk) + 1))$x)

############################################################
## Import functions that determine RF variable importance ##
############################################################

source("RandomForestFunc.R")

#################
## Run RF code ##
#################

set.seed(224672)
res_full <- run_rf()
save(res_full, file = "saves/2019-12-16-rf-full.RData")
load("saves/2019-12-16-rf-full.RData")

plot(res_full)

set.seed(645372)
res_lab_tissue <- run_rf(var_groups = c("lab", "tissue"))
save(res_lab_tissue, file = "saves/2019-12-16-rf-lab-tissue.RData")
load("saves/2019-12-16-rf-lab-tissue.RData")

plot(res_lab_tissue)

set.seed(364627)
res_celltype_lab <- run_rf(var_groups = c("celltype", "lab"))
save(res_celltype_lab, file = "saves/2019-12-16-rf-celltype-lab.RData")
load("saves/2019-12-16-rf-celltype-lab.RData")

plot(res_celltype_lab)

set.seed(292837)
res_sc_or_bulk_tissue <- run_rf(var_groups = c("sc_or_bulk", "tissue"))
save(
  res_sc_or_bulk_tissue,
  file = "saves/2019-12-16-rf-sc_or_bulk-tissue.RData")
load("saves/2019-12-16-rf-sc_or_bulk-tissue.RData")

plot(res_sc_or_bulk_tissue)

set.seed(985734)
res_subcelltype_tissue <- run_rf(var_groups = c("subcelltype", "tissue"))
save(
  res_subcelltype_tissue,
  file = "saves/2019-12-16-rf-subcelltype-tissue.RData")
load("saves/2019-12-16-rf-subcelltype-tissue.RData")

plot(res_subcelltype_tissue)

set.seed(738837)
res_bulk_lab_tissue <- run_rf(
  var_groups = c("lab", "tissue"),
  rows = which(sc_or_bulk == 1))
save(
  res_bulk_lab_tissue,
  file = "saves/2019-12-18-rf-bulk-lab-tissue.RData")
load("saves/2019-12-18-rf-bulk-lab-tissue.RData")

plot(res_bulk_lab_tissue, suffix = "bulk")

set.seed(128562)
res_sc_lab_tissue <- run_rf(
  var_groups = c("lab", "tissue"),
  rows = which(sc_or_bulk == 0))
save(
  res_sc_lab_tissue,
  file = "saves/2019-12-18-rf-sc-lab-tissue.RData")
load("saves/2019-12-18-rf-sc-lab-tissue.RData")

plot(res_sc_lab_tissue, suffix = "sc")

set.seed(87463)
res_lab_4_6_celltype <- run_rf(
  var_groups = c("celltype", "lab"),
  rows = which(lab %in% c(4, 6)),
  remove_const_cols = TRUE)
save(
  res_lab_4_6_celltype,
  file = "saves/2019-12-18-rf-sc-lab-4-6-celltype.RData")
load("saves/2019-12-18-rf-sc-lab-4-6-celltype.RData")

plot(res_lab_4_6_celltype, suffix = "lab-4-6")

set.seed(74658)
res_celltype_tissue <- run_rf(var_groups = c("celltype", "tissue"))
save(res_celltype_tissue, file = "saves/2019-12-18-rf-celltype-tissue.RData")
load("saves/2019-12-18-rf-celltype-tissue.RData")

plot(res_celltype_tissue)

set.seed(263856)
res_celltype_tissue_lab_5_bulk <- run_rf(
  var_groups = c("celltype", "tissue"),
  rows = which(lab == 5 & sc_or_bulk == 1),
  remove_const_cols = TRUE)
save(
  res_celltype_tissue_lab_5_bulk,
  file = "saves/2019-12-18-rf-bulk-lab-5-celltype-tissue.RData")
load("saves/2019-12-18-rf-bulk-lab-5-celltype-tissue.RData")

plot(res_celltype_tissue_lab_5_bulk, suffix = "bulk-lab-5")

set.seed(573626)
res_celltype_tissue_lab_5_bulk_downsampled <- run_rf(
  var_groups = c("celltype", "tissue"),
  rows = which(lab == 5 & sc_or_bulk == 1),
  remove_const_cols = TRUE,
  stratified_bootstrap = tissue[which(lab == 5 & sc_or_bulk == 1)],
  bootstrap_samples = 6)
save(
  res_celltype_tissue_lab_5_bulk_downsampled,
  file = "saves/2019-12-18-rf-bulk-lab-5-celltype-tissue-downsampled.RData")
load("saves/2019-12-18-rf-bulk-lab-5-celltype-tissue-downsampled.RData")

plot(
  res_celltype_tissue_lab_5_bulk_downsampled,
  suffix = "bulk-lab-5-downsampled")

set.seed(484738)
res_subcelltype_tissue_lab_5_bulk_downsampled <- run_rf(
  var_groups = c("subcelltype", "tissue"),
  rows = which(lab == 5 & sc_or_bulk == 1),
  remove_const_cols = TRUE,
  stratified_bootstrap = tissue[which(lab == 5 & sc_or_bulk == 1)],
  bootstrap_samples = 6)
save(
  res_subcelltype_tissue_lab_5_bulk_downsampled,
  file = "saves/2019-12-19-rf-bulk-lab-5-subcelltype-tissue-downsampled.RData")
load("saves/2019-12-19-rf-bulk-lab-5-subcelltype-tissue-downsampled.RData")

plot(
  res_celltype_tissue_lab_5_bulk_downsampled,
  suffix = "bulk-lab-5-downsampled")

set.seed(736453)
res_tissue_lab_5_bulk_downsampled <- run_rf(
  var_groups = c("subcelltype"),
  rows = which(lab == 5 & sc_or_bulk == 1),
  remove_const_cols = TRUE,
  stratified_bootstrap = tissue[which(lab == 5 & sc_or_bulk == 1)],
  bootstrap_samples = 6)
save(
  res_tissue_lab_5_bulk_downsampled,
  file = "saves/2019-12-19-rf-bulk-lab-5-tissue-downsampled.RData")
load("saves/2019-12-19-rf-bulk-lab-5-tissue-downsampled.RData")

plot(
  res_celltype_tissue_lab_5_bulk_downsampled,
  suffix = "bulk-lab-5-downsampled")

set.seed(17274)
res_sc_or_bulk <- run_rf(var_groups = c("sc_or_bulk"))
save(res_sc_or_bulk, file = "saves/2019-12-19-rf-sc_or_bulk.RData")
load("saves/2019-12-19-rf-sc_or_bulk.RData")

plot(res_sc_or_bulk)

####################
## Model checking ##
####################

## For these I would have to calculate variable importance per response variable
## again. Maybe implement at a later point
# tibble(
#   variable = as.factor(colnames(X)),
#   var_imp = oob_var_imp[,1]) %>%
#   ggplot() +
#   geom_bar(aes(x = variable, y = var_imp), stat = "identity") +
#   labs(title = "VI for PC1", x = NULL, y = "Variable Importance") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ggsave("plots/2019-11-26-PC1-var-imp.png", device = "png",
#   width = 15, height = 5)

# tibble(
#   variable = as.factor(colnames(X)),
#   var_imp = vmp[,2]) %>%
#   ggplot() +
#   geom_bar(aes(x = variable, y = var_imp), stat = "identity") +
#   labs(title = "VI for PC2", x = NULL, y = "Variable Importance") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ggsave("plots/2019-11-26-PC2-var-imp.png", device = "png",
#   width = 15, height = 5)

### Residuals

# Calculate predictions on whole data set
Y_pred <- lapply(res, function(r) {
  Reduce(`+`, lapply(r$forest, function(tree) {
    do.call(cbind, lapply(predict(tree, as.data.frame(X))$regrOutput, function(x) x$predicted))
  }))
})

residuals <- mapply(function(Y_, pred) {
  as.vector(unname(Y_) - pred)
}, Y, Y_pred, SIMPLIFY = FALSE)

ggpubr::ggarrange(plotlist = mapply(function(pred, res, type) {
  tibble(
    idx = 1L:(nrow(pred) * ncol(pred)),
    res = res) %>%
    ggplot() +
    geom_vline(
      aes(xintercept = x),
      data = tibble(x = seq(
        0.5, nrow(pred) * ncol(pred) + 0.5,
        by = nrow(pred))), size = 0.2, alpha = 0.5) +
    geom_point(aes(x = idx, y = res), size = 0.6) +
    labs(
      title = sprintf("Residuals by variable in Y (%s normalised)", type),
      x = "Index", y = "Residual")
}, Y_pred, residuals, c("TPM", "TMM", "Quantile"), SIMPLIFY = FALSE),
  nrow = 3, ncol = 1)

ggsave("plots/2019-12-06-residuals-by-var.png", device = "png",
  width = 8, height = 10)


ggpubr::ggarrange(plotlist = mapply(function(pred, res, type) {
  tibble(
    idx = 1L:(nrow(pred) * ncol(pred)),
    res = res[sample(nrow(pred) * ncol(pred))]) %>%
    ggplot() +
    geom_vline(
      aes(xintercept = x),
      data = tibble(x = seq(
        0.5, nrow(pred) * ncol(pred) + 0.5,
        by = nrow(pred))), size = 0.2, alpha = 0.5) +
    geom_point(aes(x = idx, y = res), size = 0.6) +
    labs(
      title = sprintf("Residuals in randomized order (%s normalised)", type),
      x = "Index", y = "Residual")
}, Y_pred, residuals, c("TPM", "TMM", "Quantile"), SIMPLIFY = FALSE),
  nrow = 3, ncol = 1)

ggsave("plots/2019-12-06-residuals-randomized.png", device = "png",
  width = 8, height = 10)

ggpubr::ggarrange(plotlist = mapply(function(pred, res, type) {
  tibble(
    idx = 1L:(nrow(pred) * ncol(pred)),
    res = res[sample(nrow(pred) * ncol(pred))]) %>%
    ggplot() +
    stat_qq(aes(sample = res), size = 1) +
    labs(
      title =
        sprintf("Empirical vs normal distribution (%s normalised)", type))
}, Y_pred, residuals, c("TPM", "TMM", "Quantile"), SIMPLIFY = FALSE),
  nrow = 1, ncol = 3)

ggsave("plots/2019-12-06-residuals-qqplot.png", device = "png",
  width = 18, height = 5)
