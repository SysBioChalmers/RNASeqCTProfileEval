#This script collects both bulk and single-cell data and puts them into a common matrix.
#The bulk is previously collected and processed using GenerateBulkDataMatrix.R
#The single-cell is generated with CreateSingleCellProfiles.m


# Required packages: (need to be installed)
# - edgeR (called from TMMNorm)
# - preprocessCore (for quantile normalisation)
# - readxl (for reading the design matrix)
# - sva (for ComBat batch correction)

###########
## Setup ##
###########

# Adjust to folder where data is stored
data_folder = "./data/"

# Load the MakeTPM helper function while avoiding to clutter
# the global environment
func_env <- new.env(parent = globalenv())
source("FigureHelpFunc.R", local = func_env)
MakeTPM = func_env$MakeTPM
TMMNorm = func_env$TMMNorm

###############
## Load data ##
###############
# First read TPM and count matrix
tpm_bulk = read.table(
    file = paste0(data_folder, "tpmMatrix.txt"), 
    header = TRUE, sep = "\t")

# Don't read the TMM from file here. Instead, calculate it by scaling the TPM
#  to the same library size as the counts. Then, rescale all samples the same,
# to an average of 10^6.
# bulk_tmm = read.table(
#   file = paste0(data_folder, "tmmMatrix.txt"),
#   header = TRUE, sep = "\t")
counts_bulk = read.table(
    file = paste0(data_folder, "countsMatrix.txt"), 
    header = TRUE, sep = "\t")

total_counts = colSums(counts_bulk);
pseudo_counts_bulk = t(t(tpm_bulk / 10 ^ 6) * total_counts);
# colSums(pseudo_counts_bulk) / total_counts # test

# Read the single-cell pooled samples
uc_sc = read.table(
    file = paste0(data_folder, "scProfiles.txt"),
    header = TRUE, sep = "\t", row.names = 1)

# TPM normalisation
tpm_sc = MakeTPM(uc_sc)
# colSums(tpm_sc) # test

# Merge two data frames by ID
pseudo_counts_sc_and_bulk = merge(pseudo_counts_bulk, uc_sc, by = "row.names")
# Save row names and remove corresponding column
row.names(pseudo_counts_sc_and_bulk) = pseudo_counts_sc_and_bulk$Row.names
pseudo_counts_sc_and_bulk = pseudo_counts_sc_and_bulk[,-1]

# TPM normalisation
tpm_sc_and_bulk = MakeTPM(pseudo_counts_sc_and_bulk)
# colSums(tpm_sc_and_bulk) #test

tpm_sc_and_bulk_non_filtered = tpm_sc_and_bulk
# Filter all lowly expressed genes
sel = rowMeans(tpm_sc_and_bulk) > 1
tpm_sc_and_bulk = tpm_sc_and_bulk[sel,]

tmm_sc_and_bulk_non_filtered = TMMNorm(pseudo_counts_sc_and_bulk)
tmm_sc_and_bulk = tmm_sc_and_bulk_non_filtered[sel,]

quantile_sc_and_bulk = preprocessCore::normalize.quantiles(
    as.matrix(tpm_sc_and_bulk))
quantile_sc_and_bulk_non_filtered = preprocessCore::normalize.quantiles(
    as.matrix(tpm_sc_and_bulk_non_filtered))

# Read the design matrix
dm_in <- readxl::read_xlsx(
    paste0(data_folder, "DesignMatrix.xlsx"),
    sheet = "DesignMatrix", range = "B3:DC10",
    col_names = c("variable", 1:105),
    col_types = c("text", rep.int("numeric", 105)))
  
dm_mat <- t(as.matrix(dm_in[,-1]))
mode(dm_mat) <- "integer"
colnames(dm_mat) <- dm_in$variable
dm <- as.data.frame(dm_mat)

# Do batch correction with combat
log_fold_change_tmm_sc_and_bulk = as.matrix(log2(tmm_sc_and_bulk + 0.05))
celltype = dm$`Cell type B=1` - 1 # 0 or 1 depending on B cell or T cell
batch = dm$`Lab`

mm_combat = model.matrix(
    ~ 1 + celltype, data = as.data.frame(celltype))
combat_corrected_data = sva::ComBat(
    dat = log_fold_change_tmm_sc_and_bulk, 
    batch = batch, mod = mm_combat, 
    par.prior = TRUE, prior.plots = FALSE)

# combat_pcs <- prcomp(t(combat_corrected_data))$x[,1:2]
# plot(combat_pcs) # test, looks reasonable

# Transform back
bc_sc_and_bulk = 2 ^ combat_corrected_data - 0.05
bc_sc_and_bulk[bc_sc_and_bulk < 0] = 0

# Write results to files:
saveRDS(tpm_sc_and_bulk, paste0(data_folder, "tpmScAndBulk.RDS"))
saveRDS(tmm_sc_and_bulk, paste0(data_folder, "tmmScAndBulk.RDS"))
saveRDS(quantile_sc_and_bulk, paste0(data_folder, "quantileScAndBulk.RDS"))
saveRDS(bc_sc_and_bulk, paste0(data_folder, "bcScAndBulk.RDS"))

saveRDS(
    tpm_sc_and_bulk_non_filtered, 
    paste0(data_folder, "tpmScAndBulkNonFilt.RDS"))
saveRDS(
    tmm_sc_and_bulk_non_filtered, 
    paste0(data_folder, "tmmScAndBulkNonFilt.RDS"))
saveRDS(
    quantile_sc_and_bulk_non_filtered, 
    paste0(data_folder, "quantileScAndBulkNonFilt.RDS"))

saveRDS(dm$`Bulk=1`, paste0(data_folder, "scOrBulk.RDS"))
saveRDS(dm$`Cell type B=1`, paste0(data_folder, "cellTypes.RDS"))
saveRDS(dm$`Lab`, paste0(data_folder, "labs.RDS"))
saveRDS(dm$`Sub Cell type`, paste0(data_folder, "subCellTypes.RDS"))
saveRDS(dm$`Tissue`, paste0(data_folder, "tissues.RDS"))
saveRDS(dm$`Same Individual`, paste0(data_folder, "individual.RDS"))
saveRDS(dm$`Technical Replicates`, paste0(data_folder, "techRepl.RDS"))
