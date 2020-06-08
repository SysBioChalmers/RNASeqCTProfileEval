if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("variancePartition")
browseVignettes("variancePartition")

###############
# The vignette
###############

# load library
library('variancePartition')
# load simulated data:
# geneExpr: matrix of gene expression values
# info: information/metadata about each sample
data(varPartData)
# Specify variables to consider
# Age is continuous so model it as a fixed effect
# Individual and Tissue are both categorical,
# so model them as random effects
# Note the syntax used to specify random effects
form <- ~ Age + (1|Individual) + (1|Tissue) + (1|Batch)
# Fit model and extract results
# 1) fit linear mixed model on gene expression
# If categorical variables are specified,
# a linear mixed model is used
# If all variables are modeled as fixed effects,
# a linear model is used
# each entry in results is a regression model fit on a single gene
# 2) extract variance fractions from each model fit

# for each gene, returns fraction of variation attributable
# to each variable
# Interpretation: the variance explained by each variables
# after correcting for all other variables
# Note that geneExpr can either be a matrix,
# and EList output by voom() in the limma package,
# or an ExpressionSet
varPart <- fitExtractVarPartModel( geneExpr, form, info )
# sort variables (i.e. columns) by median fraction
# of variance explained
vp <- sortCols( varPart )
# Figure 1a
# Bar plot of variance fractions for the first 10 genes
plotPercentBars( vp[1:10,] )
#
# Figure 1b
# violin plot of contribution of each variable to total variance
plotVarPart( vp )


###############
# Our project
###############



