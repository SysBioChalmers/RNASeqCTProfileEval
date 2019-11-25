library("ggplot2")
library("Seurat")

dataFolder = "C:/Work/R/RNASeqCTProfileEval/"
source(paste0(dataFolder, "FigureHelpFunc.R"))


#############################
# load data
#############################
library(dplyr)
library(tibble)

#First read TPM, tmm and count matrix
bulk_tpm = read.table(file=paste0(dataFolder, "/tpmMatrix.txt"), header=T, sep="\t")
#Don't read the tmm from file here. Instead, calculate it by scaling the TPM to the same library size
#as the counts. Then, rescale all samples the same, to an average of 10^6.

#bulk_tmm = read.table(file=paste0(dataFolder, "/tmmMatrix.txt"), header=T, sep="\t")
bulk_counts = read.table(file=paste0(dataFolder, "/countsMatrix.txt"), header=T, sep="\t")

totCounts = colSums(bulk_counts);
bulkPseudoCounts = t(t(bulk_tpm / 10^6) * totCounts);

#colSums(bulkPseudoCounts) / totCounts # test

#read the single-cell pooled samples
sc_uc = read.table(file=paste0(dataFolder, "/scProfiles.txt"), header=T, sep="\t", row.names=1)
sc_tpm = MakeTPM(sc_uc);
#colSums(sc_tpm) #test

# merge two data frames by ID
pseudoCountsScAndBulk <- merge(bulkPseudoCounts, sc_uc, by="row.names")
row.names(pseudoCountsScAndBulk) = pseudoCountsScAndBulk$Row.names
pseudoCountsScAndBulk = pseudoCountsScAndBulk[,-1]

tpmScAndBulk = MakeTPM(pseudoCountsScAndBulk)
colSums(tpmScAndBulk) #test

tpmScAndBulkNonFilt = tpmScAndBulk
#filter all lowly expressed genes
sel = rowMeans(tpmScAndBulk) > 1
tpmScAndBulk = tpmScAndBulk[sel,]


#BiocManager::install("preprocessCore")
library("edgeR")

tmmScAndBulkNonFilt = TMMNorm(pseudoCountsScAndBulk)

tmmScAndBulk = tmmScAndBulkNonFilt[sel,]


library("preprocessCore")
quantileScAndBulk = normalize.quantiles(as.matrix(tpmScAndBulk))
quantileScAndBulkNonFilt = normalize.quantiles(as.matrix(tpmScAndBulkNonFilt))



#read the design matrix
#install.packages("xlsx")
library("xlsx")
desm <- read.xlsx(paste0(dataFolder, "/DesignMatrix.xlsx"), sheetName = "DesignMatrix")
cellTypes = as.numeric(desm[7, 2:106])
labs = as.numeric(desm[2, 2:106])
subCellTypes = as.numeric(desm[6, 2:106])
tissues = as.numeric(desm[3, 2:106])
individual = as.numeric(desm[5, 2:106])
techRepl = as.numeric(desm[8, 2:106])


#do batch correction with combat
library(sva)
logtmms = as.matrix(log2(tmmScAndBulk + 0.05))
ctVar = cellTypes - 1;#0 or 1 depending on b cell or t cell
batch = labs;

modcombat = model.matrix(~1 + ctVar, data=as.data.frame(ctVar))
combat_edata = ComBat(dat=logtmms, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

#transform back
bcScAndBulk = combat_edata^2 - 0.05
bcScAndBulk[bcScAndBulk<0] = 0

#plotPCA(as.matrix(bcScAndBulk), col=as.numeric(cellTypes)) #for test, looks reasonable

#Now write to files:


saveRDS(tpmScAndBulk, paste0(dataFolder, "tpmScAndBulk.RDS"))
saveRDS(tmmScAndBulk, paste0(dataFolder, "tmmScAndBulk.RDS"))
saveRDS(quantileScAndBulk, paste0(dataFolder, "quantileScAndBulk.RDS"))
saveRDS(bcScAndBulk, paste0(dataFolder, "bcScAndBulk.RDS"))

saveRDS(tpmScAndBulkNonFilt, paste0(dataFolder, "tpmScAndBulkNonFilt.RDS"))
saveRDS(tmmScAndBulkNonFilt, paste0(dataFolder, "tmmScAndBulkNonFilt.RDS"))
saveRDS(bcScAndBulk, paste0(dataFolder, "quantileScAndBulkNonFilt.RDS"))

saveRDS(cellTypes, paste0(dataFolder, "cellTypes.RDS"))
saveRDS(labs, paste0(dataFolder, "labs.RDS"))
saveRDS(subCellTypes, paste0(dataFolder, "subCellTypes.RDS"))
saveRDS(tissues, paste0(dataFolder, "tissues.RDS"))
saveRDS(individual, paste0(dataFolder, "individual.RDS"))
saveRDS(techRepl, paste0(dataFolder, "techRepl.RDS"))



