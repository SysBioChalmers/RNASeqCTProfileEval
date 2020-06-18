#Generates Fig 4.

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
  tpm = prcomp(log2(t(tpm_sc_and_bulk) + 0.05))$x,
  tmm = prcomp(log2(t(tmm_sc_and_bulk) + 0.05))$x,
  quantile = prcomp(log2(t(quantile_sc_and_bulk) + 0.05))$x)

############################################################
## Import functions that determine RF variable importance ##
############################################################

source("RandomForestFunc.R")

###########################
## Chosen for manuscript ##
###########################

### Comment out to avoid re-running the random forest estimation

set.seed(263856)
res_celltype_tissue_lab_5_bulk <- run_rf(
  var_groups = c("celltype", "tissue"),
  rows = which(lab == 5 & sc_or_bulk == 1),
  remove_const_cols = TRUE)
save(
  res_celltype_tissue_lab_5_bulk,
  file = "saves/2019-12-19-rf-bulk-lab-5-celltype-tissue.RData")
# load("saves/2019-12-19-rf-bulk-lab-5-celltype-tissue.RData")

plot(res_celltype_tissue_lab_5_bulk, suffix = "bulk-lab-5")

set.seed(923963)
res_subcelltype_tissue_lab_5_bulk_only_b <- run_rf(
  var_groups = c("subcelltype", "tissue"),
  rows = which(lab == 5 & sc_or_bulk == 1 & celltype == 1),
  remove_const_cols = TRUE)
save(
  res_subcelltype_tissue_lab_5_bulk_only_b,
  file = paste0(
    "saves/2019-12-19-","rf-bulk-lab-5-only-b-subcelltype-tissue.RData"))
# load(paste0(
#   "saves/2019-12-19-","rf-bulk-lab-5-only-b-subcelltype-tissue.RData"))

plot(
  res_subcelltype_tissue_lab_5_bulk_only_b,
  suffix = "bulk-lab-5-only-b")

set.seed(262648)
res_subcelltype_tissue_lab_5_bulk_only_t <- run_rf(
  var_groups = c("subcelltype", "tissue"),
  rows = which(lab == 5 & sc_or_bulk == 1 & celltype == 2),
  remove_const_cols = TRUE)
save(
  res_subcelltype_tissue_lab_5_bulk_only_t,
  file = paste0(
    "saves/2019-12-19-","rf-bulk-lab-5-only-t-subcelltype-tissue.RData"))
# load(paste0(
#   "saves/2019-12-19-","rf-bulk-lab-5-only-t-subcelltype-tissue.RData"))

plot(
  res_subcelltype_tissue_lab_5_bulk_only_t,
  suffix = "bulk-lab-5-only-t")

set.seed(7563537)
res_celltype_lab_bulk <- run_rf(
  var_groups = c("celltype", "lab"),
  rows = which(tissue == 1 & subcelltype %in% c(0, 1) & sc_or_bulk == 1),
  remove_const_cols = TRUE)
save(
  res_celltype_lab_bulk,
  file = "saves/2019-12-19-","rf-bulk-celltype-lab.RData")
# load("saves/2019-12-19-rf-bulk-celltype-lab.RData")

plot(res_celltype_lab_bulk, suffix = "bulk")

set.seed(1245252)
res_celltype_lab_sc <- run_rf(
  var_groups = c("celltype", "lab"),
  rows = which(tissue == 1 & subcelltype %in% c(0, 1) & sc_or_bulk == 0),
  remove_const_cols = TRUE)
save(
  res_celltype_lab_sc,
  file = "saves/2019-12-19-rf-sc-celltype-lab.RData")
# load("saves/2019-12-19-rf-sc-celltype-lab.RData")

plot(res_celltype_lab_sc, suffix = "sc")

set.seed(947463)
res_sc_or_bulk_celltype_lab_5_6 <- run_rf(
  var_groups = c("sc_or_bulk", "celltype"),
  rows = which(tissue == 1 & lab %in% c(5, 6)),
  remove_const_cols = TRUE)
save(
  res_sc_or_bulk_celltype_lab_5_6,
  file = "saves/2019-12-19-rf-lab-5-6-sc_or_bulk-celltype.RData")
# load("saves/2019-12-19-rf-lab-5-6-sc_or_bulk-celltype.RData")

plot(res_sc_or_bulk_celltype_lab_5_6, suffix = "lab-5-6")

#######################
## Create final plot ##
#######################

 load("saves/2019-12-19-rf-bulk-lab-5-celltype-tissue.RData")
 load("saves/2019-12-19-rf-bulk-lab-5-only-b-subcelltype-tissue.RData")
 load("saves/2019-12-19-rf-bulk-lab-5-only-t-subcelltype-tissue.RData")
 load("saves/2019-12-19-rf-bulk-celltype-lab.RData")
 load("saves/2019-12-19-rf-sc-celltype-lab.RData")
 load("saves/2019-12-19-rf-lab-5-6-sc_or_bulk-celltype.RData")

create_subplot <- function(res, var_groups, title, ylab = FALSE) {
  oob_mean_var_groups_imp <- lapply(
    res, function(r) colMeans(r$oob_var_groups_imp))

  p <- tibble(
      variable = factor(rep(var_groups, times = 3), levels = var_groups),
      var_imp = do.call(c,
        lapply(oob_mean_var_groups_imp, function(vi) vi / sum(vi))),
      norm = factor(
        rep(c("TPM", "TMM", "Quantile"), each = length(var_groups)),
        levels = c("TPM", "TMM", "Quantile"))) %>%
      ggplot() +
      geom_bar(
        aes(x = variable, y = var_imp, fill = norm),
        stat = "identity",
        position = "dodge") +
      labs(x = NULL, y = NULL, title = title) +
      scale_y_continuous(limits = c(0, 1)) +
      scale_fill_manual("Normalisation", values = cbPalette) +
      theme_minimal() +
      theme(plot.title = element_text(size = 10, face = "bold"))

  if (ylab) {
    p <- p +
      labs(x = NULL, y = "Variable Importance (in %)") +
      theme(plot.margin = unit(c(0.2, 0, 0.2, 0.55), "cm"),
        axis.title = element_text(size = 8))
  } else {
    p <- p + theme(plot.margin = unit(c(0.2, 0, 0.2, 1), "cm"))
  }

  # Return final plot
  p
}

p1 <- create_subplot(
  res_sc_or_bulk_celltype_lab_5_6,
  c("Mix of sc and bulk", "Cell type"),
  "Mix of Single-Cell and Bulk vs Cell Type")

p2 <- create_subplot(
  res_celltype_tissue_lab_5_bulk,
  c("Cell type", "Tissue"),
  "Cell Type vs Tissue")

p3 <- create_subplot(
  res_subcelltype_tissue_lab_5_bulk_only_b,
  c("Cell subtype", "Tissue"),
  "Cell Subtype vs Tissue, B cells",
  ylab = TRUE)

p4 <- create_subplot(
  res_subcelltype_tissue_lab_5_bulk_only_t,
  c("Cell subtype", "Tissue"),
  "Cell Subtype vs Tissue, T cells")

p5 <- create_subplot(
  res_celltype_lab_bulk,
  c("Cell type", "Lab"),
  "Cell Type vs Lab, Bulk")

p6 <- create_subplot(
  res_celltype_lab_sc,
  c("Cell type", "Lab"),
  "Cell Type vs Lab, Single-Cell")

fig4 = ggpubr::ggarrange(
  p1, p2, p3, p4, p5, p6,
  ncol = 2, nrow = 3,
  common.legend = TRUE, legend = "bottom",
  labels = c("A", "B", "C", "D", "E", "F"),
  label.x = 0.03)

library(ggpubr)

#check that the title is shown on the graph, it sometimes randomly disappears. 
annotate_figure(fig4,
                top = text_grob("Estimating Variation Factors Using Random Forest Regression", face = "bold", size = 14))


ggsave("plots/rf-var-group-importance.png",
  width = 20.32, height = 18, unit = "cm", dpi = 100)

ggsave("plots/rf-var-group-importance.svg",
  width = 20.32, height = 18, unit = "cm", dpi = 100)

####################
## Model checking ##
####################

## Residuals

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
