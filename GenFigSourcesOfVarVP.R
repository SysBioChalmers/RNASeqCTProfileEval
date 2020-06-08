#Here, the variancePartition package is used!


#Run Generate MixedDataMatrices before running this file!
library("ggplot2")
library("Seurat")
library("ggpubr")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("variancePartition")
library('variancePartition')

fig_path = "Z:/projects/Cell type profiles/figures/"


dataFolder = "C:/Work/R/RNASeqCTProfileEval/"
source(paste0(dataFolder, "FigureHelpFunc.R"))


#############################
# load data
#############################
library(dplyr)
library(tibble)

tmmScAndBulk = readRDS(paste0(dataFolder, "data/tmmScAndBulk.RDS"))
quantileScAndBulk = readRDS(paste0(dataFolder, "data/quantileScAndBulk.RDS"))

tmmScAndBulkFilt = tmmScAndBulk[rowMeans(tmmScAndBulk) >=10,]

#log transform
geneExprLog = log2(tmmScAndBulk + 1)
geneExprFilt = log2(tmmScAndBulkFilt + 1)
geneExprLogQ = log2(quantileScAndBulk + 1)


cellTypes = readRDS(paste0(dataFolder, "data/cellTypes.RDS"))
cellTypesF = factor(cellTypes, c(1,2), c("B cell", "T cell"))
labs = readRDS(paste0(dataFolder, "data/labs.RDS"))
labsF = as.factor(labs)
subCellTypes = readRDS(paste0(dataFolder, "data/subCellTypes.RDS"))
subCellTypesF = as.factor(subCellTypes)
tissues = readRDS(paste0(dataFolder, "data/tissues.RDS"))
tissuesF = as.factor(tissues)
individual = readRDS(paste0(dataFolder, "data/individual.RDS"))
individualF = as.factor(individual)
techRepl = readRDS(paste0(dataFolder, "data/techRepl.RDS"))
techReplF = as.factor(techRepl)
scOrBulk = as.factor(readRDS(paste0(dataFolder, "data/scOrBulk.RDS")))

meta = data.frame(cellType = cellTypesF, lab=labsF, subCellType = subCellTypesF, tissue = tissuesF, scOrBulk=scOrBulk, individual = individualF, techRepl= techReplF, stringsAsFactors = F)

#Filter out so we only have one sample from the same individual
sel = c(-3,-4,-6,-7,-8,-10,-11,-12,-14,-15,-17,-18,-20,-21,-23,-24)
diffIndExpr = geneExprLog[, sel]
diffIndMeta = meta[sel,]
metaBulkDI = diffIndMeta[diffIndMeta$scOrBulk == 1,]
exprBulkDI = diffIndExpr[, diffIndMeta$scOrBulk == 1]
exprBulkQDI = geneExprLogQ[, sel][, diffIndMeta$scOrBulk == 1]

sel2 = c(sel, -104, -105)
exprBulkDIDSSC = geneExprLog[,sel2]
metaBulkDIDSSC = meta[sel2,]


metaBulk = meta[meta$scOrBulk == 1,]
exprBulk = geneExprLog[, meta$scOrBulk == 1]
exprBulkQ = geneExprLogQ[, meta$scOrBulk == 1]
exprBulkFilt = geneExprFilt[, meta$scOrBulk == 1]

exprSC = geneExprLog[, meta$scOrBulk == 0]
exprSCQ = geneExprLogQ[, meta$scOrBulk == 0]
metaSC = meta[meta$scOrBulk == 0,]

exprNoSmartSeq = geneExprLog[,1:103]
metaNoSmartSeq = meta[1:103,]

#Drop-based single cell (i.e. skip smart-seq)
exprDSSC = geneExprLog[, meta$scOrBulk == 0 & meta$lab != 10]#some rows are 0, which doesn't work
exprDSSC = exprDSSC[rowSums(exprDSSC) != 0,]
exprDSSCQ = geneExprLogQ[, meta$scOrBulk == 0  & meta$lab != 10]
exprDSSCQ = exprDSSCQ[rowSums(exprDSSCQ) != 0,]
metaDSSC = meta[meta$scOrBulk == 0 & meta$lab != 10,]


##########################################
#Create plots for lab, cell subtype and tissue for different gene sets


form <- ~ (1|subCellType) + (1|tissue) + (1|lab)

varPartDgs <- fitExtractVarPartModel(exprBulkDI, form, metaBulkDI)
colnames(varPartDgs) = c("Lab", "Cell subtype", "Tissue", "Residuals")

#All genes
pAllGenes = plotVarPart(varPartDgs, main="All Genes", ylab = "Variance expl. (%)")
pAllGenes
dim(exprBulkDI)#12072

#Housekeeping genes:
hkGenes = read.table(paste0(dataFolder, "data/HK_genes.txt"), stringsAsFactors = F)$V1
varPartHKG = varPartDgs[row.names(varPartDgs) %in% hkGenes, ]
dim(varPartHKG)#3393 genes
pHKGenes = plotVarPart(varPartHKG, main = "HK Genes", ylab = "Variance expl. (%)")
pHKGenes

#Lm22 genes:
lm22 = read.table(paste0(dataFolder, "data/LM22.txt"), header = T, sep="\t", stringsAsFactors = F)
lm22genes = lm22$Gene.symbol
varPartLm22 = varPartDgs[row.names(varPartDgs) %in% lm22genes, ]
dim(varPartLm22)#395 genes

pLm22Genes = plotVarPart(varPartLm22, main = "LM22 Genes", ylab = "Variance expl. (%)")
pLm22Genes

#now, select the genes in lm22 that differ much between B and T cells

#the columns in lm22 that are B and T cells
bSel = 2:4 #the cell types that are some kind of B cells
tSel = 5:11 #same for T cells

rmb = rowMeans(lm22[,bSel])
rmt = rowMeans(lm22[,tSel])

abslfc = abs(log2(rmb/rmt))
lm22highDiffGenes = lm22genes[abslfc > 1] #So, all genes where the lfc between B and T is at least 1

varPartLm22HighDiff = varPartDgs[row.names(varPartDgs) %in% lm22highDiffGenes, ]
dim(varPartLm22HighDiff)#274 genes

pLm22GenesLFC = plotVarPart(varPartLm22HighDiff, main = "LM22 LFC(B,T) > 1", ylab = "Variance expl. (%)")
pLm22GenesLFC

#Plot when using cell type instead of cell subtype
form <- ~ (1|cellType) + (1|tissue) + (1|lab)

varPartDgsCT <- fitExtractVarPartModel(exprBulkDI, form, metaBulkDI)
varPartDgsCTSrt = varPartDgsCT[, c(2,1,3,4)]
colnames(varPartDgsCTSrt) = c("Lab", "Cell type", "Tissue", "Residuals")
varPartDgsCTSrtFilt = varPartDgsCTSrt[row.names(varPartDgsCTSrt) %in% lm22highDiffGenes, ]
pCT = plotVarPart(varPartDgsCTSrtFilt, main = "Cell Type (B/T)", ylab = "Variance expl. (%)")
pCT

#And, finally, plot the explained variance as a function of gene expression:
#test to plot per gene expression
meanGeneExpression = rowMeans(exprBulkDI)

sortx = sort(meanGeneExpression, index.return=T)

ds = cbind(x=sortx$x, varPartDgs[sortx$ix,])

loess_fit <- loess(Lab~x, ds, span = 0.3)
predLab = predict(loess_fit)
loess_fit <- loess(`Cell subtype`~x, ds, span = 0.3)
predSubCellType = predict(loess_fit)
loess_fit <- loess(Tissue~x, ds, span = 0.3)
predTissue = predict(loess_fit)
loess_fit <- loess(Residuals~x, ds, span = 0.3)
predResiduals = predict(loess_fit)

dsPlot = data.frame(x=ds$x, Lab=predLab, "Cell Subtype"=predSubCellType, Tissue=predTissue, Residuals = predResiduals)
#dsPlot = data.frame(x=ds$x, Lab=predLab,Tissue=predTissue, Residuals = predResiduals)

library(reshape2)
dsPlot.m <- melt(dsPlot, id='x')
levels(dsPlot.m$variable) = c("Lab", "Cell subtype", "Tissue", "Residuals")
colnames(dsPlot.m)[2] = "Factor"

pPerGEX = ggplot(data=dsPlot.m, aes(x=x, y=value, group=Factor, colour = Factor)) +
  geom_line() +
  ggtitle("Expl. Var. per Gene Expr.") +
  xlab(expression(Log[2]*"(pseudo-TPM)")) + 
  ylab("Variation expl.") +
  scale_colour_manual(values=c(ggColorHue(3), "grey65")) +
  #theme(legend.position= "bottom", legend.direction = "vertical", legend.box = "horizontal", legend.title = element_blank()) +
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5)) +
  theme(text 		= element_text(colour="black"), 
        axis.text 	= element_text(colour="black"),
        legend.text = element_text(colour="black")) +
  #theme(legend.position= "bottom", legend.direction = "vertical", legend.box = "horizontal", legend.title = element_blank()) +
  theme(legend.position= "bottom", legend.direction = "vertical", legend.title = element_blank(), legend.margin=margin(t = -0.6, unit='cm')) +
  theme(plot.title = element_text(lineheight=.8, face="bold")) +
  guides(colour=guide_legend(ncol=2, nrow=2, byrow=F))

pPerGEX


figGeneSetsComb = ggarrange(pAllGenes,pHKGenes,pLm22Genes,pLm22GenesLFC,pCT, pPerGEX, nrow=3, ncol=2, labels=c("A","B","C","D","E","F"))

figGeneSets = annotate_figure(figGeneSetsComb,
                       top = text_grob("Explained Variance for Bulk Samples", face = "bold", size = 14))
figGeneSets

ggsave(
  paste0(fig_path, "Fig3.png"),
  plot = figGeneSets, device = "png",
  width = 6, height = 9, dpi = 300)


################################
## Single-cell variation figure
################################
colors = c(ggColorHue(3), "grey85")
formDSSC = ~ (1|cellType) + (1|lab) + (1|tissue)
varPartDSSC = fitExtractVarPartModel(exprDSSC, formDSSC, metaDSSC)
#resort and rename columns
varPartDSSCSrt = varPartDSSC[,c(2,1,3,4)]
colnames(varPartDSSCSrt) = c("Lab", "Cell type", "Tissue", "Residuals")

#all genes
pDSSCAllGenes = plotVarPart(varPartDSSCSrt, main="All Genes", ylab = "Variance expl. (%)", col = colors)
pDSSCAllGenes

#LM22, LFC > 1
varPartDSSCLm22HighDiff = varPartDSSCSrt[row.names(varPartDSSCSrt) %in% lm22highDiffGenes, ]
pDSSCLm22GenesLFC = plotVarPart(varPartDSSCLm22HighDiff, main = "LM22 LFC(B,T) > 1", ylab = "Variance expl. (%)", col = colors)
pDSSCLm22GenesLFC

#mixed, LM22, LFC > 1
varBulkDIDSSC = fitExtractVarPartModel(exprBulkDIDSSC, formDSSC, metaBulkDIDSSC)
#resort and rename columns
varBulkDIDSSCSrt = varBulkDIDSSC[,c(2,1,3,4)]
colnames(varBulkDIDSSCSrt) = c("Lab", "Cell type", "Tissue", "Residuals")
varPartBulkDIDSSCLm22HighDiff = varBulkDIDSSCSrt[row.names(varBulkDIDSSCSrt) %in% lm22highDiffGenes, ]
pBulkDIDSSCLab = plotVarPart(varPartBulkDIDSSCLm22HighDiff, main = "SC/Bulk Mix - Lab", ylab = "Variance expl. (%)", col = colors)
pBulkDIDSSCLab

formDSSC2 = ~ (1|cellType) + (1|scOrBulk) + (1|tissue)
varBulkDIDSSC2 = fitExtractVarPartModel(exprBulkDIDSSC, formDSSC2, metaBulkDIDSSC)
#resort and rename columns
varBulkDIDSSCSrt2 = varBulkDIDSSC2[,c(2,1,3,4)]
colnames(varBulkDIDSSCSrt2) = c("SC or bulk", "Cell type", "Tissue", "Residuals")
varPartBulkDIDSSCLm22HighDiff2 = varBulkDIDSSCSrt2[row.names(varBulkDIDSSCSrt2) %in% lm22highDiffGenes, ]
pBulkDIDSSCBvsSC = plotVarPart(varPartBulkDIDSSCLm22HighDiff2, main = "SC/Bulk Mix - SC or Bulk", ylab = "Variance expl. (%)", col = colors)
pBulkDIDSSCBvsSC


figSCComb = ggarrange(pDSSCAllGenes,pDSSCLm22GenesLFC,pBulkDIDSSCLab,pBulkDIDSSCBvsSC, nrow=2, ncol=2, labels=c("A","B","C","D"))

figSC = annotate_figure(figSCComb,
                              top = text_grob("Explained Variance for Single-Cell Samples", face = "bold", size = 14))
figSC

ggsave(
  paste0(fig_path, "Fig4.png"),
  plot = figSC, device = "png",
  width = 6, height = 6, dpi = 300)



######################################
# Supplementary figure
######################################

# TMM
form = ~ (1|cellType) + (1|lab) + (1|tissue) 
varPartAllTmm = fitExtractVarPartModel(geneExprLog, form, meta)
#resort and rename columns
varPartAllTmmSrt = varPartAllTmm[,c(2,1,3,4)]
colnames(varPartAllTmmSrt) = c("Lab", "Cell type", "Tissue", "Residuals")
pAllTMM = plotVarPart(varPartAllTmmSrt, main="TMM Normalization", ylab = "Variance expl. (%)")
pAllTMM

#Quantile
form = ~ (1|cellType) + (1|lab) + (1|tissue) 
varPartAllQ = fitExtractVarPartModel(geneExprLogQ, form, meta)
#resort and rename columns
varPartAllQSrt = varPartAllQ[,c(2,1,3,4)]
colnames(varPartAllQSrt) = c("Lab", "Cell type", "Tissue", "Residuals")
pAllQ = plotVarPart(varPartAllQSrt, main="Quantile Normalization", ylab = "Variance expl. (%)")
pAllQ

#Sc with smart-seq2:
colors = c(ggColorHue(3)[-3], "grey85")
formSC = ~ (1|cellType) + (1|lab) 
varPartSCAll = fitExtractVarPartModel(exprSC, formSC, metaSC)
#resort and rename columns
varPartSCAllSrt = varPartSCAll[,c(2,1,3)]
colnames(varPartSCAllSrt) = c("Lab", "Cell type", "Residuals")
pSCAll = plotVarPart(varPartSCAllSrt, main="Single-Cell Incl. Smart-Seq2", ylab = "Variance expl. (%)", col = colors)
pSCAll

#...and without
varPartDSSC2 = fitExtractVarPartModel(exprDSSC, formSC, metaDSSC)
#resort and rename columns
varPartDSSC2Srt = varPartDSSC2[,c(2,1,3)]
colnames(varPartDSSC2Srt) = c("Lab", "Cell type", "Residuals")
pDSSC2 = plotVarPart(varPartDSSC2Srt, main="Single-Cell Excl. Smart-Seq2", ylab = "Variance expl. (%)", col = colors)
pDSSC2

#Individual in lab 4 (not technical replicates)
exprInd = geneExprLog[, labs==4]
totExprInd = rowSums(exprInd)
exprInd = exprInd[totExprInd != 0,]


metaInd = meta[labs==4,]

formInd = ~ (1|cellType) + (1|individual) 
varPartInd = fitExtractVarPartModel(exprInd, formInd, metaInd)

#resort and rename columns
varPartIndSrt = varPartInd
colnames(varPartIndSrt) = c("Cell type", "Individual", "Residuals")
colorsInd = c(ggColorHue(3)[-3], "grey85")
colorsInd[1] = colorsInd[2]
colorsInd[2] = rgb(1,0,1)
pInd = plotVarPart(varPartIndSrt, main="Individual", ylab = "Variance expl. (%)", col = colorsInd)
pInd


#technical replicates in lab 3
exprTech = geneExprLog[, labs==3]
totExprTech = rowSums(exprTech)
exprTech = exprTech[totExprTech != 0,]

metaTech = meta[labs==3,]

formTech = ~ (1|techRepl) 
varPartTech = fitExtractVarPartModel(exprTech, formTech, metaTech)

#resort and rename columns
varPartTechSrt = varPartTech
colnames(varPartTechSrt) = c("Techn. Repl.", "Residuals")
colorsTech = c(ggColorHue(3)[c(-2, -3)], "grey85")
colorsTech[1] = rgb(0,1,1)
pTech = plotVarPart(varPartTechSrt, main="Technical Replicates", ylab = "Variance expl. (%)", col = colorsTech)
pTech

#Metabolic genes (not so interesting, skip):
metGenes = read.table(paste0(dataFolder, "data/metabolic_genes.txt"), stringsAsFactors = F)$V1
varPartMet = varPartDgs[row.names(varPartDgs) %in% metGenes, ]
dim(varPartMet)#2679 genes
pMetGenes = plotVarPart(varPartMet, main = "Metabolic Genes", ylab = "Variance expl. (%)")
pMetGenes

figS1Comb = ggarrange(pAllTMM,pAllQ,pSCAll,pDSSC2,pInd,pTech, nrow=3, ncol=2, labels=c("A","B","C","D","E","F"))

figS1 = annotate_figure(figS1Comb,
                        top = text_grob("Suppl. Figures for Explained Variance", face = "bold", size = 14))
figS1

ggsave(
  paste0(fig_path, "FigS1.png"),
  plot = figS1, device = "png",
  width = 6, height = 9, dpi = 300)


#Check for outliers
#######################

#bulk first
numBulk = dim(exprBulk)[2]
varPartRes = vector(mode = "list", length = numBulk)
form <- ~ (1|subCellType) + (1|tissue) + (1|lab)
for (i in 1:numBulk) {
  varPartRes[[i]] = fitExtractVarPartModel(exprBulk[,-i], form, metaBulk[-i,])
}

varPartRes
saveRDS(varPartRes, "outlierTestBulk.RDS")

varPartBulkFull = fitExtractVarPartModel(exprBulk, form, metaBulk)

labTot = mean(varPartBulkFull$lab)
sctTot = mean(varPartBulkFull$subCellType)
tisTot = mean(varPartBulkFull$tissue)
resTot = mean(varPartBulkFull$Residuals)
labOutl = rep(0,length(varPartRes));
sctOutl = rep(0,length(varPartRes));
tisOutl = rep(0,length(varPartRes));
resOutl = rep(0,length(varPartRes));
for (i in 1:length(varPartRes)) {
  labOutl[i] = mean(varPartRes[[i]]$lab)
  sctOutl[i] = mean(varPartRes[[i]]$subCellType)
  tisOutl[i] = mean(varPartRes[[i]]$tissue)
  resOutl[i] = mean(varPartRes[[i]]$Residuals)
}
labDiffs = -labOutl + labTot
labDiffsAbs = abs(labDiffs)
sctDiffsAbs = abs(-sctOutl + sctTot)
tisDiffsAbs = abs(-tisOutl + tisTot)
resDiffsAbs = abs(-resOutl + resTot)

sort(labDiffsAbs)
sort(sctDiffsAbs)
sort(tisDiffsAbs)
sort(tisDiffsAbs)

#The maximum change is less than 2 percentage units for any of the factors

#now, for the single-cell samples

numSc = dim(exprDSSC)[2]
varPartResSC = vector(mode = "list", length = numSc)
form <- ~ (1|cellType) + (1|tissue) + (1|lab)
for (i in c(1:23, 25:numSc)) {

  toRem = -i
  #there is a problem with i = 23, since there is only one sample left with tissue = 5.
  #so, for this case, remove both 23 and 24
  if (i == 23 ) {
    toRem = c(-i, -i-1)
  }
  exprDSSCTmp = exprDSSC[,toRem]
  exprDSSCTmp = exprDSSCTmp[rowSums(exprDSSCTmp) != 0,]
  varPartResSC[[i]] = fitExtractVarPartModel(exprDSSCTmp, form, metaDSSC[toRem,])
  if (i == 23) {
    varPartResSC[[i+1]] = varPartResSC[[i]];
  }
}

varPartResSC
saveRDS(varPartResSC, "outlierTestSC.RDS")

varPartSCFull = fitExtractVarPartModel(exprDSSC, form, metaDSSC)

labTot = mean(varPartSCFull$lab)
ctTot = mean(varPartSCFull$cellType)
tisTot = mean(varPartSCFull$tissue)
resTot = mean(varPartSCFull$Residuals)
labOutl = rep(0,length(varPartResSC));
ctOutl = rep(0,length(varPartResSC));
tisOutl = rep(0,length(varPartResSC));
resOutl = rep(0,length(varPartResSC));
for (i in 1:length(varPartResSC)) {
  labOutl[i] = mean(varPartResSC[[i]]$lab)
  ctOutl[i] = mean(varPartResSC[[i]]$cellType)
  tisOutl[i] = mean(varPartResSC[[i]]$tissue)
  resOutl[i] = mean(varPartResSC[[i]]$Residuals)
}
labDiffs = -labOutl + labTot
labDiffsAbs = abs(labDiffs)
ctDiffsAbs = abs(-ctOutl + ctTot)
tisDiffsAbs = abs(-tisOutl + tisTot)
resDiffsAbs = abs(-resOutl + resTot)

sort(labDiffsAbs)
sort(sctDiffsAbs)
sort(tisDiffsAbs)
sort(tisDiffsAbs)

#The maximum change is less than 2 percentage units for any of the factors here as well!
