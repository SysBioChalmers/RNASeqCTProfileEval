#also try RUV:
#BiocManager::install("RUVSeq")
#library("RUVSeq")

#read the design matrix
#install.packages("xlsx")
library("xlsx")
dm <- read.xlsx("C:/Work/R/RNASeqCTProfileEval/DesignMatrix.xlsx", sheetName = "DesignMatrix")
cellTypes = as.numeric(dm[9, 2:75])
labs = as.numeric(dm[2, 2:75])
subCellTypes = as.numeric(dm[7, 2:75])
tissues = as.numeric(dm[3, 2:75])


#get housekeeping genes
pathHkGenes = "C:/Work/MatlabCode/components/SCLib/ImportableData/HK_genes.txt"
hkGenes = read.table(pathHkGenes)
hkGenes = hkGenes[,1]


library(RUVSeq)
filter <- apply(counts, 1, function(x) length(x[x>5])>=2)
filtered <- counts[filter,]

#genes <- rownames(filtered)[grep("^ENS", rownames(filtered))]
#spikes <- rownames(filtered)[grep("^ERCC", rownames(filtered))]

x <- as.factor(cellTypes)

library(RColorBrewer)
colors <- brewer.pal(3, "Set2")
plotRLE(counts, outline=FALSE, ylim=c(-4, 4), col=as.numeric(cellTypes))
plotPCA(counts, col=as.numeric(cellTypes), cex=1.2)

uc = betweenLaneNormalization(filtered, which="upper")
plotRLE(uc, outline=FALSE, ylim=c(-4, 4), col=as.numeric(cellTypes))
plotPCA(uc, col=as.numeric(cellTypes), cex=1.2)

plotRLE(tmms, outline=FALSE, ylim=c(-4, 4), col=as.numeric(cellTypes))
plotPCA(tmms, col=as.numeric(cellTypes), cex=1.2)


plotRLE(tpms, outline=FALSE, ylim=c(-4, 4), col=as.numeric(cellTypes))


#it seems that upper quartile normalization does not do a very good job, compared to TMM?


RUVInput = newSeqExpressionSet(as.matrix(uc),
                               phenoData = data.frame(x, row.names=colnames(uc)))

#RUVInput = newSeqExpressionSet(as.matrix(uc),
#                               phenoData = data.frame(x, row.names=colnames(tmms)))


hkGenes2 = intersect(hkGenes, row.names(uc))

RUVed = RUVg(RUVInput, hkGenes2, k=1)
plotPCA(RUVed, col=cellTypes, cex=1.2)
plotPCA(RUVed, col=labs, cex=1.2)

#test combat

#BiocManager::install("sva")
#library("sva")
#library(bladderbatch)
#data(bladderdata)
#library(pamr)
#library(limma)

#first, filter lowly expressed genes (based on counts)
filter <- apply(counts, 1, function(x) length(x[x>5])>=2)
filttmms = tmms[filter,]

#create log-transformed data from TMM normalization
#logtmms = unlist(as.data.frame(log2(filttmms + 1)))
logtmms = as.matrix(log2(filttmms + 1))
ctVar = cellTypes - 1;#0 or 1 depending on b cell or t cell
batch = labs;

modcombat = model.matrix(~1 + ctVar, data=as.data.frame(ctVar))

#modcombat2 = model.matrix(~1, data=pheno)
#edata = exprs(bladderEset)

#batch2 = pheno$batch
#combat_edata2 = ComBat(dat=edata, batch=batch2, mod=modcombat2, par.prior=TRUE, prior.plots=FALSE)


combat_edata = ComBat(dat=logtmms, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

#try pca on corrected data
plotPCA((2^combat_edata)-1, col=as.numeric(cellTypes), cex=1.2)
plotPCA((2^combat_edata)-1, col=as.numeric(subCellTypes), cex=1.2)
plotPCA((2^combat_edata)-1, col=as.numeric(tissues), cex=1.2)
plotPCA((2^combat_edata)-1, col=as.numeric(labs), cex=1.2)

plotPCA(tmms, col=as.numeric(labs), cex=1.2)



if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
BiocManager::install(c("Biobase","sva","bladderbatch","snpStats"))
library(devtools)
library(Biobase)
library(sva)
library(bladderbatch)
library(snpStats)
data(bladderdata)
