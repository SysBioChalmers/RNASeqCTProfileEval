library(ggplot2)
library(dplyr)
library(randomForestSRC)
# The following additional packages need to be installed:
# - stringr
# - scales
# - ggpubr

data_folder <- "./data/"

cbPalette <- c(
	"#999999", "#E69F00", "#56B4E9", "#009E73",
	"#F0E442", "#0072B2", "#D55E00", "#CC79A7") # colour-blind friendly palette

### Consider random forests
#
# Random forests can naturally deal with regression of approximating y by
# a functional relationship f(x) without requiring that f is linear.
# Multivariate random forests are of particular interest here since they can
#
# For now choose variables
# - Lab
# - Tissue
# - Bulk or SC
# - Cell Type
# - Sub Cell Type

##################
## Convert data ##
##################

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

var_groups <- c("sc_or_bulk", "celltype", "subcelltype", "lab", "tissue")

# Load gene datasets (samples in columns)
tpm_sc_and_bulk <- readRDS(paste0(data_folder, "tpmScAndBulk.RDS"))
tmm_sc_and_bulk <- readRDS(paste0(data_folder, "tmmScAndBulk.RDS"))
quantile_sc_and_bulk <- readRDS(paste0(data_folder, "quantileScAndBulk.RDS"))

# Exclude technical replicates
tech_reps <- c(3, 4, 6, 7, 11)
X <- X[-tech_reps,]
tpm_sc_and_bulk <- tpm_sc_and_bulk[,-tech_reps]
tmm_sc_and_bulk <- tmm_sc_and_bulk[,-tech_reps]
quantile_sc_and_bulk <- quantile_sc_and_bulk[,-tech_reps]

tibble(
  variable = factor(colnames(X), levels = colnames(X)),
  frequency = colSums(X)) %>%
  ggplot() +
  geom_bar(aes(x = variable, y = frequency), stat = "identity") +
  labs(title = "Frequency of covariates", x = NULL, y = "Frequency") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("plots/2019-12-06-covariate-frequency.png", device = "png",
  width = 15, height = 5)

##### Work in full principal component space, i.e. no information loss, just
##### working in 100 dimensional space instead of ~ 15000
Y <- list(
  tpm = prcomp(log2(t(tpm_sc_and_bulk) + 1))$x,
  tmm = prcomp(log2(t(tmm_sc_and_bulk) + 1))$x,
  quantile = prcomp(log2(t(quantile_sc_and_bulk) + 1))$x)

###########################################################
## Homebrewed solution using the randomForestSRC package ##
###########################################################
set.seed(224672)

# Settings
n_tree <- 300
n_feature <- floor(dim(X)[2] / 3) # Use standard recommendation for number
                                  # of features for regression (#features / 3)
min_leaf <- 5 # Standard setting

# Number of times the permutations for variable importance are repeated and averaged over
# This is inspired by
# Ishwaran H (2007) Variable importance in binary regression trees and forests.
# Electronic Journal of Statistics 1:519-537 DOI 10.1214/07-EJS039
n_perms <- 5

res <- lapply(c("tpm", "tmm", "quantile"), function(type) {

n_samples <- nrow(Y[[type]])
n_outputs <- ncol(Y[[type]])

fml <- as.formula(paste0(
  "Multivar(", paste(sprintf("PC%d", 1:n_outputs), collapse = ", "), ") ~ ."))

# Create bootstrap samples to be used to build the regression trees
bs_sets_res <- bootstrap::bootstrap(1:n_samples, n_tree, identity)
bs_sets <- bs_sets_res$thetastar

n_var <- ncol(X)

oob_error <- rep.int(0, n_tree) # Overall OOB error (calculated at each iter)
oob_error_tree <- rep.int(0, n_tree) # Tree-specific OOB error
oob_pred <- list() # Out-of-bag predictions for each tree

# OOB-based variable importance
# Since variables were one-hot encoded whole groups of pseudo-variables (all
# belonging to one of the original variables) will be permuted.
#
# Computational idea: Permute variable group j and re-calculate
# OOB error of the tree; compare that error to the tree-specific OOB.
# Repeat procedure multiple times.
#
# Large differences in outcome (positive and negative) stand for variable
# groups that were important in that specific tree.
# Should be aggregated over trees in the end.
oob_var_groups_imp <- matrix(0, nrow = n_tree, ncol = length(var_groups))

bs_test_set <- lapply(
  1:n_tree, function(i) (1:n_samples)[-sort(unique(bs_sets[,i]))])

forest <- list()

for (i in 1:n_tree) {
  cat(sprintf(
    paste0("%s - Build tree #%", ceiling(log10(n_tree)), "d"),
    as.POSIXct(Sys.time()), i))
  X_tr <- X[bs_sets[,i],, drop = FALSE]
  Y_tr <- Y[[type]][bs_sets[,i],, drop = FALSE]

  X_te <- X[bs_test_set[[i]],, drop = FALSE]
  Y_te <- Y[[type]][bs_test_set[[i]],, drop = FALSE]

  tree <- rfsrc.cart(
    fml, data = as.data.frame(cbind(Y_tr, X_tr)),
    mtry = n_feature, nodesize = min_leaf)

  forest[[i]] <- tree

  # Out-of-bag predictions
  oob_pred[[i]] <- do.call(
    cbind,
    lapply(predict(tree, as.data.frame(X_te))$regrOutput, function(x) x$predicted))

  # Out-of-bag error for this specific tree
  oob_error_tree[i] <-  mean((Y_te - oob_pred[[i]]) ^ 2)

  # Variable group permutations
  for (j in 1L:length(var_groups)) {
    for (l in 1:n_perms) {
      X_te_perm <- X_te
      var_idx <- stringr::str_starts(colnames(X_te_perm), var_groups[j])
      row_idx <- sample(1:nrow(X_te_perm))
      X_te_perm[,var_idx] <- X_te_perm[row_idx, var_idx]
      oob_var_groups_imp[i,j] <- (oob_var_groups_imp[i,j] + mean(
        (Y_te - do.call(
            cbind,
            lapply
            (predict(tree, as.data.frame(X_te_perm))$regrOutput,
            function(x) x$predicted))) ^ 2) -
        oob_error_tree[i])
    }
  }
  oob_var_groups_imp[i,] <- oob_var_groups_imp[i,] / n_perms

  # Overall (mean) OOB error
  oob_error[i] <- mean(do.call(c, lapply(1:n_samples, function(k) {
    # Determine the trees where sample k was not used during training
    tree_mask <- colSums(bs_sets[, 1:i, drop = FALSE] == k) == 0
    # If there were any, then calculate their average prediction and compare
    # to the observed sample. Otherwise simply return NULL
    if (any(tree_mask)) {
      preds <- do.call(rbind, lapply(
        which(tree_mask), function(j) {
          row_idx <- which(bs_test_set[[j]] == k)
          oob_pred[[j]][row_idx,]
        }))
      mean((Y[[type]][k,] - colMeans(preds)) ^ 2)
    } else {
      NULL
    }
  })))

  cat(sprintf(" OOB Error: %1.6e\n", oob_error[i]))
}

# Return summary about the random forest
list(
  forest = forest,
  oob_error = oob_error,
  oob_var_groups_imp = oob_var_groups_imp,
  bootstrap = list(
    bs_sets = bs_sets,
    bs_test_set = bs_test_set),
  oob = list(
    oob_error_tree = oob_error_tree,
    oob_pred = oob_pred))
}) # End of lapply over dataset type

save(res, file = "saves/2019-12-06-rf-full-run.RData")

## Analysis
ggpubr::ggarrange(plotlist = mapply(function(r, type) {
  tibble(err = r$oob_error, ind = 1L:length(r$oob_error)) %>%
    ggplot() +
    geom_line(aes(x = ind, y = err)) +
    labs(
      title = sprintf("RF on %s normalised data", type),
      x = "Iteration",
      y = "OOB Error")
}, res, c("TPM", "TMM", "Quantile"), SIMPLIFY = FALSE),
  nrow = 3, ncol = 1)

ggsave("plots/2019-12-06-oob-error-convergence-plot.png",
  width = 5, height = 6)

oob_mean_var_groups_imp <- lapply(
  res, function(r) colMeans(r$oob_var_groups_imp))
for(i in 1L:length(oob_mean_var_groups_imp)) {
  names(oob_mean_var_groups_imp[[i]]) <- var_groups
}

tibble(
  variable = factor(
    rep(
      c("bulk vs sc", "cell type", "sub cell type", "lab", "tissue"),
      times = 3),
    levels = c("bulk vs sc", "cell type", "sub cell type", "lab", "tissue")),
  var_imp = do.call(c, oob_mean_var_groups_imp),
  norm = factor(
    rep(c("TPM", "TMM", "Quantile"), each = 5),
    levels = c("TPM", "TMM", "Quantile"))) %>%
  ggplot() +
  geom_bar(
    aes(x = variable, y = var_imp, fill = norm),
    stat = "identity",
    position = "dodge") +
  labs(x = NULL, y = "Variable Importance") +
  scale_fill_manual("Normalisation", values = cbPalette)

ggsave("plots/2019-12-06-var-imp-per-group.png", device = "png",
  width = 8, height = 2.5)

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
