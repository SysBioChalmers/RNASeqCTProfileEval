library(tidyverse)
library(scales)
library(randomForestSRC)

source("R/data-prep-final.R")

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
#
# and perform one-hot encoding on the categorical variables
# `Lab`, `Tissue` and `Sub Cell Type`

##################
## Convert data ##
##################

one_hot_enc <- function(data, old_var, new_var) {
    enc <- do.call(
        cbind, lapply(sort(unique(data[old_var]) %>% as_vector()),
                      function(i) as.integer(data[old_var] == i)))
    colnames(enc) <- sprintf(
      "%s%d", new_var,
      sort(unique(data[old_var]) %>% as_vector()))

    # Return encoded vector
    enc
}

bulk_enc <- one_hot_enc(dm, "Bulk=1", "bulk")
celltype_enc <- one_hot_enc(dm, "Cell type B=1", "celltype")
subcelltype_enc <- one_hot_enc(dm, "Sub Cell type", "subcelltype")
lab_enc <- one_hot_enc(dm, "Lab", "lab")
tissue_enc <- one_hot_enc(dm, "Tissue", "tissue")

# Encoded design matrix
X <- bind_cols(
    as_tibble(bulk_enc), as_tibble(celltype_enc),
    as_tibble(subcelltype_enc),
    as_tibble(lab_enc), as_tibble(tissue_enc)) %>%
    as.matrix()

# Exclude technical replicates
X <- X[-c(3, 4, 6, 7, 11),]
bulk_sc_tpm_mat <- bulk_sc_tpm_mat[-c(3, 4, 6, 7, 11),]

# Split up those samples that belong to sub-cell type 1 (mix) by whether they are
# B or T cells
idx <- rowSums(X[,c("celltype1", "subcelltype1")]) == 2
X[idx, "subcelltype1"] <- 0
nms <- colnames(X)
nms <- c(nms[1:4], "subcelltype0", nms[5:30])
X <- cbind(X[,1:4], as.numeric(idx), X[,5:30])
colnames(X) <- nms

tibble(
  variable = factor(colnames(X), levels = colnames(X)),
  frequency = colSums(X)) %>%
  ggplot() +
  geom_bar(aes(x = variable, y = frequency), stat = "identity") +
  labs(title = "Frequency of covariates", x = NULL, y = "Frequency") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("plots/2019-11-27-covariate-frequency.png", device = "png",
  width = 6, height = 5)

##### Work in full principal component space, i.e. no information loss, just
##### working in 100 dimensional space instead of ~ 15000
Y <- prcomp(log2(bulk_sc_tpm_mat + 1))$x


###########################################################
## Homebrewed solution using the randomForestSRC package ##
###########################################################
set.seed(224672)

fml <- as.formula(paste0(
  "Multivar(", paste(sprintf("PC%d", 1:100), collapse = ", "), ") ~ ."))

n_samples <- nrow(Y)
n_outputs <- ncol(Y)

# Somewhat following the help for `variable_importance_measure` in
# package `MultivariateRandomForest`
n_tree <- 150
n_feature <- floor(dim(X)[2] / 3) # Use standard recommendation for number
                                  # of features for regression (#features / 3)
min_leaf <- 5 # Standard setting

# Create bootstrap samples to be used to build the regression trees
bs_sets_res <- bootstrap::bootstrap(1:n_samples, n_tree, identity)
bs_sets <- bs_sets_res$thetastar

n_var <- ncol(X)

oob_error <- rep.int(0, n_tree) # Overall OOB error (calculated at each iter)
oob_error_tree <- rep.int(0, n_tree) # Tree-specific OOB error
oob_pred <- list() # Out-of-bag predictions for each tree

# OOB-based variable importance, i.e. permute variable j and re-calculate
# OOB error of the tree; compare that error to the tree-specific OOB.
# Large differences in outcome (positive and negative) stand for variables
# that were important in that specific tree.
# Should be aggregated over trees in the end.
oob_var_imp <- matrix(0, nrow = n_tree, ncol = n_var)
oob_var_imp_joint <- list()

# Number of times the permutations for variable importance are repeated and averaged over
# This is inspired by
# Ishwaran H (2007) Variable importance in binary regression trees and forests.
# Electronic Journal of Statistics 1:519-537 DOI 10.1214/07-EJS039
n_perms <- 5

bs_test_set <- lapply(
  1:n_tree, function(i) (1:n_samples)[-sort(unique(bs_sets[,i]))])

forest <- list()

for (i in 1:n_tree) {
  cat(sprintf(
    paste0("%s - Build tree #%", ceiling(log10(n_tree)), "d"),
    as.POSIXct(Sys.time()), i))
  tree <- NULL
  X_tr <- X[bs_sets[,i],, drop = FALSE]
  Y_tr <- Y[bs_sets[,i],, drop = FALSE]

  X_te <- X[bs_test_set[[i]],, drop = FALSE]
  Y_te <- Y[bs_test_set[[i]],, drop = FALSE]

  tree <- rfsrc.cart(
    fml, data = as_tibble(cbind(Y_tr, X_tr)), 
    mtry = n_feature, nodesize = min_leaf)

  forest[[i]] <- tree

  # Out-of-bag predictions
  oob_pred[[i]] <- do.call(
    cbind,
    lapply(predict(tree, as_tibble(X_te))$regrOutput, function(x) x$predicted))

  # Out-of-bag error for this specific tree
  oob_error_tree[i] <-  mean((Y_te - oob_pred[[i]]) ^ 2)

  # One-variable permutation
  for (j in 1:n_var) {
    for (l in 1:n_perms) {
      X_te_perm <- X_te
      X_te_perm[,j] <- sample(X_te_perm[,j])
      oob_var_imp[i,j] <- (oob_var_imp[i,j] + mean(
        (Y_te - do.call(
            cbind,
            lapply
            (predict(tree, as_tibble(X_te_perm))$regrOutput, 
            function(x) x$predicted))) ^ 2) -
        oob_error_tree[i])
    }
  }
  oob_var_imp[i,] <- oob_var_imp[i,] / n_perms

  # Two-variable permutations
  oob_var_imp_joint[[i]] <- matrix(0, nrow = n_var, ncol = n_var)
  for (j in 1:(n_var - 1)) {
    for (k in (j + 1):n_var) {
      for (l in 1:n_perms) {
        X_te_perm <- X_te
        X_te_perm[,j] <- sample(X_te_perm[,j])
        X_te_perm[,k] <- sample(X_te_perm[,k])
        oob_var_imp_joint[[i]][j,k] <- (oob_var_imp_joint[[i]][j,k] + mean(
          (Y_te - do.call(
            cbind,
            lapply
              (predict(tree, as_tibble(X_te_perm))$regrOutput, 
              function(x) x$predicted))) ^ 2) -
          oob_error_tree[i])
      }
    }
  }
  oob_var_imp_joint[[i]] <- oob_var_imp_joint[[i]] / n_perms

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
      mean((Y[k,] - colMeans(preds)) ^ 2)
    } else {
      NULL
    }
  })))

  cat(sprintf(" OOB Error: %1.6e\n", oob_error[i]))
}

save.image("saves/2019-11-26-full-run.RData")

## Analysis
plot(oob_error, type = "l")

oob_mean_var_imp <- colMeans(oob_var_imp)
names(oob_mean_var_imp) <- colnames(X)

bulk_vs_sc_var_imp <- mean(oob_mean_var_imp[1:2])
celltype_var_imp <- mean(oob_mean_var_imp[3:4])
subcelltype_var_imp <- mean(abs(oob_mean_var_imp[5:17]))
lab_var_imp <- mean(oob_mean_var_imp[18:27])
tissue_var_imp <- mean(oob_mean_var_imp[28:31])

oob_mean_var_imp_joint <- Reduce(`+`, oob_var_imp_joint) / n_tree
for (i in 1:n_var) {
  for (j in 1:i) {
    oob_mean_var_imp_joint[i, j] <- NA
  }
}

oob_mean_var_imp_add <- matrix(NA, nrow = n_var, ncol = n_var)
for (i in 1:(n_var - 1)) {
    for (j in (i + 1):n_var) {
        oob_mean_var_imp_add[i, j] <- oob_mean_var_imp[i] + oob_mean_var_imp[j]
    }
}

oob_mean_var_imp_assoc <- oob_mean_var_imp_joint - oob_mean_var_imp_add
# oob_mean_var_imp_assoc <- (oob_mean_var_imp_assoc - min(oob_mean_var_imp_assoc, na.rm = TRUE))
# oob_mean_var_imp_assoc <- oob_mean_var_imp_assoc / max(oob_mean_var_imp_assoc, na.rm = TRUE)

tibble(
  variable = factor(
    c("bulk vs sc", "cell type", "sub cell type", "lab", "tissue"),
    levels = c("bulk vs sc", "cell type", "sub cell type", "lab", "tissue")),
  var_imp = c(bulk_vs_sc_var_imp,
              celltype_var_imp, subcelltype_var_imp,
              lab_var_imp, tissue_var_imp)) %>%
  ggplot() +
  geom_bar(aes(x = variable, y = var_imp), stat = "identity") +
  labs(x = NULL, y = "Variable Importance")

ggsave("plots/2019-11-27-variable-importance-overall.png", device = "png",
  width = 4, height = 2.5)

tibble(
  variable = factor(names(oob_mean_var_imp), levels = names(oob_mean_var_imp)),
  var_imp = oob_mean_var_imp / max(oob_mean_var_imp)) %>%
  ggplot() +
  geom_bar(aes(x = variable, y = var_imp), stat = "identity") +
  labs(
    title = "Variable importance across response variables",
    x = NULL, y = "Standardised Variable Importance") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("plots/2019-11-27-variable-importance-detail.png", device = "png",
  width = 15, height = 5)

tibble(
  value = as.vector(oob_mean_var_imp_assoc),
  x = rep(n_var - 1:n_var + 1, each = n_var),
  y = rep(1:n_var, times = n_var)) %>%
  ggplot() +
  geom_tile(aes(x = x, y = y, fill = value, colour = value), width = 1, height = 1) +
  coord_fixed() +
  scale_x_continuous(
      breaks = 1:30, labels = rev(colnames(X)[-1]),
      expand = expand_scale(add = c(0.25, 0.25))) +
  scale_y_continuous(
      breaks = 1:30, labels = colnames(X)[-31],
      expand = expand_scale(add = c(0.25, 0.25))) +
  scale_fill_gradient2(low = muted("blue"), mid = "#eeeeee", high = "red", na.value = "transparent") +
  scale_colour_continuous(low = "black", high = "black", na.value = "transparent", guide = FALSE) +
  labs(x = NULL, y = NULL) +
  theme(
    panel.background = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave("plots/2019-11-27-variable-associations.png", device = "png",
    width = 8, height = 8)

tibble(
  value = as.vector(oob_mean_var_imp_joint) / max(oob_mean_var_imp_joint, na.rm = TRUE),
  x = rep(n_var - 1:n_var + 1, each = n_var),
  y = rep(1:n_var, times = n_var)) %>%
  ggplot() +
  geom_tile(aes(x = x, y = y, fill = value, colour = value), width = 1, height = 1) +
  coord_fixed() +
  scale_x_continuous(
      breaks = 1:30, labels = rev(colnames(X)[-1]),
      expand = expand_scale(add = c(0.25, 0.25))) +
  scale_y_continuous(
      breaks = 1:30, labels = colnames(X)[-31],
      expand = expand_scale(add = c(0.25, 0.25))) +
  scale_fill_gradient2(low = muted("blue"), mid = "#eeeeee", high = "red", na.value = "transparent") +
  scale_colour_continuous(low = "black", high = "black", na.value = "transparent", guide = FALSE) +
  labs(x = NULL, y = NULL) +
  theme(
    panel.background = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave("plots/2019-11-27-variable-joint.png", device = "png",
    width = 8, height = 8)