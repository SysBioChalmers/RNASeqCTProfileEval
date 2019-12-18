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

###########################################################
## Homebrewed solution using the randomForestSRC package ##
###########################################################

#' Build a random forest with variable importance for groups of variables
#'
#' Instead of calculating variable importance per variable, it is calculated
#' on a group basis.
#'
#' The repetition of permuations was inspired by
#' Ishwaran H (2007) Variable importance in binary regression trees and forests.
#' Electronic Journal of Statistics 1:519-537 DOI 10.1214/07-EJS039
#'
#' @param X The covariates. Continuous or one-hot-encoded for categorical
#'          variables (samples x features)
#' @param Y The response values (samples x responses)
#' @param var_groups A character vector of patterns. Each pattern is taken as
#'                   the start of the variable names.
#' @param ... A list of optional parameters
#'  \item{n_tree}{Number of trees (default = 300)}
#'  \item{n_feature}{Number of features randomly picked at each splitting point
#'                   (default = # features / 3)}
#'  \item{min_leaf}{Minimum number of elements in a leaf (node without children)
#'                  (default = 5)}
#'  \item{n_perms}{Number of times the permutations for variable importance
#'                 are repeated and averaged over (default = 50)}
#'
#' @returns A \code{list} with elements
#'  \item{forest}{A \code{list} of \code{rfsrc.cart} trees}
#'  \item{oob_error}{The OOB error over iterations (i.e. additional trees)}
#'  \item{oob_var_groups_imp}{Permutation variable importance determined for
#'                            each group of variables as determined by
#'                            \code{var_groups}}
#'  \item{bootstrap}{A \code{list} containing
#'    \item{bs_training_sets}{The indices for the training sets}
#'    \item{bs_test_sets}{The corresponding test sets}}
#'  \item{oob}{A \code{list} containing
#'    \item{oob_error_tree}{The OOB error for each specific tree
#'                          (not aggregated as for \code{oob_error})}
#'    \item{oob_pred}{The OOB prediction for each tree}}
run_rf_var_group_imp <- function(X, Y, var_groups, ...) {
  add_opts <- list(...)

  n_tree <- add_opts$n_tree
  if (is.null(n_tree)) {
    n_tree <- 300
  }
  n_feature <- add_opts$n_feature
  if (is.null(n_feature)) {
    n_feature <- floor(dim(X)[2] / 3)
  }
  min_leaf <- add_opts$min_leaf
  if (is.null(min_leaf)) {
    min_leaf <- 5
  }
  n_perms <- add_opts$n_perms
  if (is.null(n_perms)) {
    n_perms <- 50
  }

  n_samples <- nrow(Y)
  n_outputs <- ncol(Y)

  fml <- as.formula(paste0(
    "Multivar(", paste(sprintf("PC%d", 1:n_outputs), collapse = ", "), ") ~ ."))

  # Create bootstrap samples to be used to build the regression trees
  bs_training_sets_res <- bootstrap::bootstrap(1:n_samples, n_tree, identity)
  bs_training_sets <- bs_training_sets_res$thetastar

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

  bs_test_sets <- lapply(
    1:n_tree, function(i) (1:n_samples)[-sort(unique(bs_training_sets[,i]))])

  forest <- list()

  for (i in 1:n_tree) {
    cat(sprintf(
      paste0("%s - Build tree #%", ceiling(log10(n_tree)), "d"),
      as.POSIXct(Sys.time()), i))
    X_tr <- X[bs_training_sets[,i],, drop = FALSE]
    Y_tr <- Y[bs_training_sets[,i],, drop = FALSE]

    X_te <- X[bs_test_sets[[i]],, drop = FALSE]
    Y_te <- Y[bs_test_sets[[i]],, drop = FALSE]

    tree <- rfsrc.cart(
      fml, data = as.data.frame(cbind(Y_tr, X_tr)),
      mtry = n_feature, nodesize = min_leaf)

    forest[[i]] <- tree

    # Out-of-bag predictions
    oob_pred[[i]] <- do.call(
      cbind,
      lapply(
        predict(tree, as.data.frame(X_te))$regrOutput,
        function(x) x$predicted))

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
      tree_mask <- colSums(bs_training_sets[, 1:i, drop = FALSE] == k) == 0
      # If there were any, then calculate their average prediction and compare
      # to the observed sample. Otherwise simply return NULL
      if (any(tree_mask)) {
        preds <- do.call(rbind, lapply(
          which(tree_mask), function(j) {
            row_idx <- which(bs_test_sets[[j]] == k)
            oob_pred[[j]][row_idx,]
          }))
        mean((Y[k,] - colMeans(preds)) ^ 2)
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
      bs_training_sets = bs_training_sets,
      bs_test_sets = bs_test_sets),
    oob = list(
      oob_error_tree = oob_error_tree,
      oob_pred = oob_pred))
}

##########################################################################
## A simplified interface for running the code and automatic evaluation ##
##########################################################################

# ... are additional parameters passed to `run_rf_var_group_imp`
run_rf <- function(
    var_groups = c("sc_or_bulk", "celltype", "subcelltype", "lab", "tissue"),
    rows = 1L:nrow(X), ...) {

  res <- structure(lapply(
    c("tpm", "tmm", "quantile"),
    function(type) run_rf_var_group_imp(
      X[rows, stringr::str_starts(
        colnames(X),
        sprintf("(%s)", paste0(var_groups, collapse = "|")))],
      Y[[type]][rows,],
      var_groups = var_groups,
      ...)
  ), class = "custom_rf")

  # Save meta-data
  attr(res, "var_groups") <- var_groups
  attr(res, "rows") <- rows

  res
}

#' Plot the outcome of our custom random forest implementation
#'
#' Generates some standard plots we are interested in.
#'
#' @param x The resulting object from \code{run_rf} of class \code{custom_rf}
#' @param path Path to the folder where the plots will be saved
#' @param suffix A user-defined suffic that can be used to mark different plots
plot.custom_rf <- function(x, path = "plots/", suffix = NULL, ...) {
  var_groups <- attr(x, "var_groups")

  # Determine the template for the file path
  file_path <- paste0(
    path, format(Sys.Date(), "%Y-%m-%d"),
    "-%s-",
    paste0(var_groups, collapse = "-"))
  if(!is.null(suffix)) {
    file_path <- paste0(file_path, "-", suffix)
  }
  file_path <- paste0(file_path, ".png")

  p1 <- ggpubr::ggarrange(plotlist = mapply(function(r, type) {
    tibble(err = r$oob_error, ind = 1L:length(r$oob_error)) %>%
      ggplot() +
      geom_line(aes(x = ind, y = err)) +
      labs(
        title = sprintf("RF on %s normalised data", type),
        x = "Iteration",
        y = "OOB Error")
  }, x, c("TPM", "TMM", "Quantile"), SIMPLIFY = FALSE),
    nrow = 3, ncol = 1)

  ggsave(
    sprintf(file_path, "oob-error-convergence-plot"),
    plot = p1, device = "png", width = 5, height = 6)

  oob_mean_var_groups_imp <- lapply(
    x, function(r) colMeans(r$oob_var_groups_imp))
  for(i in 1L:length(oob_mean_var_groups_imp)) {
    names(oob_mean_var_groups_imp[[i]]) <- var_groups
  }

  p2 <- tibble(
    variable = factor(
      rep(var_groups, times = 3),
      levels = var_groups),
    var_imp = do.call(c, oob_mean_var_groups_imp),
    norm = factor(
      rep(c("TPM", "TMM", "Quantile"), each = length(var_groups)),
      levels = c("TPM", "TMM", "Quantile"))) %>%
    ggplot() +
    geom_bar(
      aes(x = variable, y = var_imp, fill = norm),
      stat = "identity",
      position = "dodge") +
    labs(x = NULL, y = "Variable Importance") +
    scale_fill_manual("Normalisation", values = cbPalette) +
    theme_minimal()

  ggsave(
    sprintf(file_path, "var-imp-per-group-unnormalised"),
    plot = p2, device = "png", width = 8, height = 2.5)

  p3 <- tibble(
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
    labs(x = NULL, y = "Relative Variable Importance (in %)") +
    scale_fill_manual("Normalisation", values = cbPalette) +
    theme_minimal()

  ggsave(sprintf(file_path, "var-imp-per-group-normalised"),
    plot = p3, device = "png", width = 8, height = 2.5)
}

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
