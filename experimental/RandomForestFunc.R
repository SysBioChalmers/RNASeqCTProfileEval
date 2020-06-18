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
#'  \item{bootstrap_replace}{Whether or not bootstrap should sample
#'                           with replacement (default = TRUE)}
#'  \item{bootstrap_samples}{Number of samples per bootstrap sample
#'                           (default = nrow(X), i.e. as many as there are
#'                           samples in the data).}
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
run_rf_var_group_imp_pre <- function(X, Y, var_groups, ...) {
  if (nrow(X) != nrow(Y)) {
    stop("X and Y need to have the same number of samples (rows)")
  }

  n_samples <- nrow(Y)
  n_outputs <- ncol(Y)

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

  bootstrap_replace <- add_opts$bootstrap_replace
  if (is.null(bootstrap_replace)) {
    bootstrap_replace <- TRUE
  }

  bootstrap_samples <- add_opts$bootstrap_samples
  if (is.null(bootstrap_samples)) {
    bootstrap_samples <- n_samples
  }

  fml <- as.formula(paste0(
    "Multivar(", paste(sprintf("PC%d", 1:n_outputs), collapse = ", "), ") ~ ."))

  # Create bootstrap samples to be used to build the regression trees
  bs_training_sets <- lapply(1L:n_tree, function(i)
      sample(1L:n_samples, bootstrap_samples, replace = bootstrap_replace))

  bs_test_sets <- lapply(
      1L:n_tree, function(i)
      (1L:n_samples)[-sort(unique(bs_training_sets[[i]]))])

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

  # Save the regression trees here
  forest <- list()

  for (i in 1:n_tree) {
    cat(sprintf(
      paste0("%s - Build tree #%", ceiling(log10(n_tree)), "d"),
      as.POSIXct(Sys.time()), i))
    X_tr <- X[bs_training_sets[[i]],, drop = FALSE]
    Y_tr <- Y[bs_training_sets[[i]],, drop = FALSE]

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
      # Determine the trees where sample k was in the test set
      tree_mask <- sapply(1:i, function(j) any(bs_test_sets[[j]] == k))
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

run_rf_var_group_imp <- compiler::cmpfun(run_rf_var_group_imp_pre)

##########################################################################
## A simplified interface for running the code and automatic evaluation ##
##########################################################################

# ... are additional parameters passed to `run_rf_var_group_imp`
run_rf <- function(
    var_groups = c("sc_or_bulk", "celltype", "subcelltype", "lab", "tissue"),
    rows = 1L:nrow(X),
    remove_const_cols = FALSE,
    ...) {

  X_ <- X[rows, stringr::str_starts(
        colnames(X),
        sprintf("(%s)", paste0(var_groups, collapse = "|")))]
  if (remove_const_cols) {
    X_ <- X_[,!(colSums(X_) == 0 | colSums(X_) == nrow(X_))]
  }

  res <- structure(lapply(
    c("tpm", "tmm", "quantile"),
    function(type) run_rf_var_group_imp(
      X_,
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
  oob_sd_var_groups_imp <- lapply(
    x, function(r) apply(r$oob_var_groups_imp, 2, function(v) sd(v[v > 0])))

  p2 <- tibble(
    variable = factor(
      rep(var_groups, times = 3),
      levels = var_groups),
    var_imp = do.call(c, oob_mean_var_groups_imp),
    var_imp_min = do.call(c,
      mapply(
        function(mvi, sdvi) (mvi - sdvi),
        oob_mean_var_groups_imp, oob_sd_var_groups_imp,
        SIMPLIFY = FALSE)),
    var_imp_max = do.call(c,
      mapply(
        function(mvi, sdvi) (mvi + sdvi),
        oob_mean_var_groups_imp, oob_sd_var_groups_imp,
        SIMPLIFY = FALSE)),
    norm = factor(
      rep(c("TPM", "TMM", "Quantile"), each = length(var_groups)),
      levels = c("TPM", "TMM", "Quantile"))) %>%
    ggplot() +
    geom_bar(
      aes(x = variable, y = var_imp, fill = norm),
      stat = "identity",
      position = "dodge") +
    geom_errorbar(
      aes(x = variable, ymin = var_imp_min, ymax = var_imp_max, fill = norm),
      position = position_dodge(width = 0.9),
      width = 0.5) +
    labs(x = NULL, y = "Variable Importance") +
    scale_fill_manual("Normalisation", values = cbPalette) +
    theme_minimal() +
    theme(legend.position = "bottom")

  ggsave(
    sprintf(file_path, "var-imp-per-group-unnormalised"),
    plot = p2, device = "png",
    width = 2 + length(var_groups), height = 2.5)

  p3 <- tibble(
    variable = factor(rep(var_groups, times = 3), levels = var_groups),
    var_imp = do.call(c,
      lapply(oob_mean_var_groups_imp, function(vi) vi / sum(vi))),
    var_imp_min = do.call(c,
      mapply(
        function(mvi, sdvi) (mvi - sdvi) / sum(mvi),
        oob_mean_var_groups_imp, oob_sd_var_groups_imp,
        SIMPLIFY = FALSE)),
    var_imp_max = do.call(c,
      mapply(
        function(mvi, sdvi) (mvi + sdvi) / sum(mvi),
        oob_mean_var_groups_imp, oob_sd_var_groups_imp,
        SIMPLIFY = FALSE)),
    norm = factor(
      rep(c("TPM", "TMM", "Quantile"), each = length(var_groups)),
      levels = c("TPM", "TMM", "Quantile"))) %>%
    ggplot() +
    geom_bar(
      aes(x = variable, y = var_imp, fill = norm),
      stat = "identity",
      position = "dodge") +
    geom_errorbar(
      aes(x = variable, ymin = var_imp_min, ymax = var_imp_max, fill = norm),
      position = position_dodge(width = 0.9),
      width = 0.5) +
    labs(x = NULL, y = "Variable Importance (in %)") +
    scale_fill_manual("Normalisation", values = cbPalette) +
    theme_minimal() +
    theme(legend.position = "bottom")

  ggsave(sprintf(file_path, "var-imp-per-group-normalised"),
    plot = p3, device = "png",
    width = 2 + length(var_groups), height = 2.75)
}
