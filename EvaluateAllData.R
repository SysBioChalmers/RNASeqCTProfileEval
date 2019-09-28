
#First read TPM, tmm and count matrix
folder = "C:/Work/R/RNASeqCTProfileEval" # to be replaced by each user
bulk_tpm = read.table(file=paste0(folder, "/tpmMatrix.txt"), header=T, sep="\t")
bulk_tmm = read.table(file=paste0(folder, "/tmmMatrix.txt"), header=T, sep="\t")
bulk_counts = read.table(file=paste0(folder, "/countsMatrix.txt"), header=T, sep="\t")

totCounts = colSums(bulk_counts);

#read the single-cell pooled samples
sc_tpm = read.table(file=paste0(folder, "/scProfiles.txt"), header=T, sep="\t", row.names=1)

# merge two data frames by ID
scAndBulk <- merge(bulk_tpm,sc_tpm, by="row.names")
row.names(scAndBulk) = scAndBulk$Row.names
scAndBulk = scAndBulk[,-1]


#test merge
#a = c(1, 2, 3)
#b = c(1, 4, 5)
#df1 = data.frame(a,b)
#row.names(df1) = c("a","b", "c")
#df2 = df1
#row.names(df2) = c("a","b", "d")
#dfm = merge(df1,df2, by="row.names")
#row.names(dfm) = dfm$Row.names
#dfm = dfm[,-1]


#read the design matrix
#install.packages("xlsx")
library("xlsx")
desm <- read.xlsx(paste0(folder, "/DesignMatrix.xlsx"), sheetName = "DesignMatrix")
cellTypes = as.numeric(desm[9, 2:104])
labs = as.numeric(desm[2, 2:104])
subCellTypes = as.numeric(desm[7, 2:104])
tissues = as.numeric(desm[3, 2:104])

plotRLE(as.matrix(bulk_counts), outline=FALSE, ylim=c(-4, 4), col=as.numeric(cellTypes))
plotRLE(as.matrix(bulk_tpm), outline=FALSE, ylim=c(-4, 4), col=as.numeric(cellTypes))
plotRLE(as.matrix(bulk_tmm), outline=FALSE, ylim=c(-4, 4), col=as.numeric(cellTypes))





library(RUVSeq)

#scb = lapply(scAndBulk, as.numeric)
#sum(is.na(as.matrix(scb)))

#mat = as.matrix(as.numeric(scAndBulk))

plotPCA(as.matrix(scAndBulk), col=as.numeric(cellTypes), cex=1.2)
#plotPCA(as.matrix(scb), col=as.numeric(cellTypes), cex=1.2)

plotPCA(as.matrix(sc_tpm), cex=1.2)
plotPCA(as.matrix(bulk_tpm), cex=1.2)


#experiment with tcga lusc
tcga_tpm = read.table(file=paste0(folder, "/tcga_lusc.txt"), header=T, sep="\t", row.names=1)
tcga_mat = as.matrix(tcga_tpm)

tcga_filt = tcga_tpm[rowMeans(tcga_tpm) > 1,]
bulk_filt = bulk_tpm[rowMeans(bulk_tpm) > 1,]


a = log2((tcga_filt[,1] + 1)/(tcga_filt[,2] + 1))
b = log2((bulk_filt[,14] + 1)/(bulk_filt[,15] + 1))

c = c(a,b)

df = as.data.frame(c)
#ggplot(data=df, aes(df$a)) +
#  geom_histogram(binwidth = 0.05) +
#  geom_histogram(aes = b, binwidth = 0.05)

num = length(a);
sel = c(rep(0,num),rep(1,length(b)));

ggplot(data=df, aes(x=c)) +
  geom_histogram(data=subset(df,sel == 0),fill = "red", alpha = 0.2, binwidth = 0.05) +
  geom_histogram(data=subset(df,sel == 1),fill = "blue", alpha = 0.2, binwidth = 0.05)



#check if batch correction gets rid of dropout bias:
#first see how it looks without correction:

#run TMM (using functions in InvestigateSingleCellBiases)
scAndBulkTPM = MakeTPM(scAndBulk)
scAndBulkTMM = TMMNorm(scAndBulkTPM)

plotPCA(as.matrix(scAndBulkTPM), col=as.numeric(cellTypes), cex=1.2)
plotPCA(as.matrix(scAndBulkTMM), col=as.numeric(cellTypes), cex=1.2)



#filter rows with too low expression
filt = rowMeans(scAndBulkTPM) >= 0.1;
filtTMMData = scAndBulkTMM[filt,]

div = log2((filtTMMData[,75] + 0.01) / (filtTMMData[,10] + 0.01))
meanExpr = (filtTMMData[,75] + filtTMMData[,10]) / 2
PlotSingleCellBias(div,meanExpr)

library(sva)



logtmms = as.matrix(log2(filtTMMData + 1))
ctVar = cellTypes - 1;#0 or 1 depending on b cell or t cell
batch = labs;

modcombat = model.matrix(~1 + ctVar, data=as.data.frame(ctVar))



#modcombat2 = model.matrix(~1, data=pheno)
#edata = exprs(bladderEset)

#batch2 = pheno$batch
#combat_edata2 = ComBat(dat=edata, batch=batch2, mod=modcombat2, par.prior=TRUE, prior.plots=FALSE)


#apply(RNA_seq_log2_D1, 1, var)!=0

combat_edata = ComBat(dat=logtmms, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

#transform back
bcData = combat_edata^2 - 1
bcData[bcData<0] = 0

plotPCA(as.matrix(bcData), col=as.numeric(cellTypes), cex=1.2)


#TMM once more
bcNorm = TMMNorm(MakeTPM(bcData))
plotPCA(as.matrix(bcNorm), col=as.numeric(cellTypes), cex=1.2)

#no norm after:
divbc = log2((bcData[,75] + 0.01) / (bcData[,10] + 0.01))
meanExprbc = (bcData[,75] + bcData[,10]) / 2
PlotSingleCellBias(divbc,meanExprbc)
#TMM after - didn't make things better:
divbc = log2((bcNorm[,75] + 0.01) / (bcNorm[,10] + 0.01))
meanExprbc = (bcNorm[,75] + bcNorm[,10]) / 2
PlotSingleCellBias(divbc,meanExprbc)

#conclusion - 
# * ComBat makes the samples more similar
# * ComBat reduces, but does not completely fix the "dropout issue", lowly expressed genes are 
#   still more lowly expressed in single-cell data






