library("ggplot2")
library("Seurat")

dataFolder = "C:/Work/R/RNASeqCTProfileEval/"
source(paste0(dataFolder, "FigureHelpFunc.R"))


#############################
# load data
#############################

#Download instructions for HCA single-cell evaluation:
# The raw data are in GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132044.
# The processed data are in the Single Cell Portal.
# For cortex: https://portals.broadinstitute.org/single_cell/study/SCP425/single-cell-comparison-cortex-data
# For Mixture: https://portals.broadinstitute.org/single_cell/study/SCP426/single-cell-comparison-mixture-data
# For PBMC: https://portals.broadinstitute.org/single_cell/study/SCP424/single-cell-comparison-pbmc-data
# You will need cortex only. Unzip everything. Then,
# copy some files:
# genes.txt -> counts/genes.tsv, genes.txt -> UMIs/genes.tsv  - According to instructions, we should have picked genes.count.txt, but it is identical but weird it seems
# counts.reads.txt -> counts/matrix.mtx, counts.umis.txt -> UMIs/matrix.mtx
# cell.names.txt -> counts/barcodes.tsv, cell.names.txt -> UMIs/barcodes.tsv

HCASCE_counts <- Read10X(data.dir = "C:/Work/MatlabCode/components/SingleCellToolbox/ImportableData/HCA_single-cell_comparison/cortex/counts")
HCASCE_UMIs <- Read10X(data.dir = "C:/Work/MatlabCode/components/SingleCellToolbox/ImportableData/HCA_single-cell_comparison/cortex/UMIs")

HCASCE_bulk1 = read.table("C:/Work/MatlabCode/components/SingleCellToolbox/ImportableData/HCA_single-cell_comparison/cortex/bulk/cortex1.genes.results",header=T, sep="\t", row.names = 1)
HCASCE_bulk2 = read.table("C:/Work/MatlabCode/components/SingleCellToolbox/ImportableData/HCA_single-cell_comparison/cortex/bulk/cortex2.genes.results",header=T, sep="\t", row.names = 1)

HCASCE_genes = row.names(HCASCE_bulk1)
HCASCE_genes2 = row.names(HCASCE_bulk2) #they are the same

#convert gene ids
library(biomaRt)
ensembl_us_west = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="uswest.ensembl.org")
listDatasets(ensembl_us_west)
ensembl = useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")

geneConvTableM <- getBM(attributes=c('ensembl_gene_id', 'mgi_symbol'), mart = ensembl)

ind = match(HCASCE_genes,geneConvTableM$ensembl_gene_id)
newGenes = geneConvTableM$mgi_symbol[ind]


HCASCE_bulkTPM = cbind(HCASCE_bulk1$TPM,HCASCE_bulk2$TPM)
HCASCE_bulkCounts = cbind(HCASCE_bulk1$expected_count,HCASCE_bulk2$expected_count)
row.names(HCASCE_bulkTPM) = newGenes;
row.names(HCASCE_bulkCounts) = newGenes;

#remove genes that could not be converted
HCASCE_bulkTPM = HCASCE_bulkTPM[!is.na(row.names(HCASCE_bulkTPM)),]
HCASCE_bulkCounts = HCASCE_bulkCounts[!is.na(row.names(HCASCE_bulkCounts)),]

HCASCE_c1SmartSeq2 = HCASCE_counts[,grep("Cortex1_Smart_seq2", colnames(HCASCE_counts))]
HCASCE_c2SmartSeq2 = HCASCE_counts[,grep("Cortex2_Smart_seq2", colnames(HCASCE_counts))]

HCASCE_smartSeq2Mn = cbind(rowSums(as.matrix(HCASCE_c1SmartSeq2)), rowSums(as.matrix(HCASCE_c2SmartSeq2)))
colnames(HCASCE_smartSeq2Mn) = c("smartseq1","smartseq2")
rm(HCASCE_c1SmartSeq2, HCASCE_c2SmartSeq2)

#merge them - don't use merge, it is super slow
library(dplyr)
library(tibble)
HCASCE_mc = inner_join(rownames_to_column(as.data.frame(HCASCE_bulkCounts)), rownames_to_column(as.data.frame(HCASCE_smartSeq2Mn)))
row.names(HCASCE_mc) = HCASCE_mc$rowname
HCASCE_mc = HCASCE_mc[,-1]
HCASCE_mc = MakeTPM(HCASCE_mc)
HCASCE_mcn = TMMNorm(HCASCE_mc)

#now, do the same with 10x: UMIs and counts compared to bulk TPM
HCASCE_u110x = as.matrix(HCASCE_UMIs[,grep("Cortex1_10xChromium", colnames(HCASCE_UMIs))])
HCASCE_u210x = as.matrix(HCASCE_UMIs[,grep("Cortex2_10xChromium", colnames(HCASCE_UMIs))])
HCASCE_c110x = as.matrix(HCASCE_counts[,grep("Cortex1_10xChromium", colnames(HCASCE_counts))])
HCASCE_c210x = as.matrix(HCASCE_counts[,grep("Cortex2_10xChromium", colnames(HCASCE_counts))])
HCASCE_u10xMn = cbind(rowSums(HCASCE_u110x), rowSums(HCASCE_u210x), rowSums(HCASCE_c110x), rowSums(HCASCE_c210x))
colnames(HCASCE_u10xMn) = c("UMI1","UMI2", "count1", "count2")
rm(HCASCE_u110x, HCASCE_u210x, HCASCE_c110x, HCASCE_c210x)


HCASCE_m10x = inner_join(rownames_to_column(as.data.frame(HCASCE_bulkTPM)), rownames_to_column(as.data.frame(HCASCE_u10xMn)))
row.names(HCASCE_m10x) = HCASCE_m10x$rowname
HCASCE_m10x = HCASCE_m10x[,-1]
HCASCE_m10xun = HCASCE_m10x
HCASCE_m10x = MakeTPM(HCASCE_m10x)
HCASCE_m10xn = TMMNorm(HCASCE_m10x)




#First read TPM, tmm and count matrix
folder = "C:/Work/R/RNASeqCTProfileEval" # to be replaced by each user
bulk_tpm = read.table(file=paste0(folder, "/tpmMatrix.txt"), header=T, sep="\t")
#Don't read the tmm from file here. Instead, calculate it by scaling the TPM to the same library size
#as the counts. Then, rescale all samples the same, to an average of 10^6.

#bulk_tmm = read.table(file=paste0(folder, "/tmmMatrix.txt"), header=T, sep="\t")
bulk_counts = read.table(file=paste0(folder, "/countsMatrix.txt"), header=T, sep="\t")

totCounts = colSums(bulk_counts);
bulkPseudoCounts = t(t(bulk_tpm / 10^6) * totCounts);

#colSums(bulkPseudoCounts) / totCounts # test

#read the single-cell pooled samples
sc_uc = read.table(file=paste0(folder, "/scProfiles.txt"), header=T, sep="\t", row.names=1)
sc_tpm = MakeTPM(sc_uc);
#colSums(sc_tpm) #test

# merge two data frames by ID
pseudoCountsScAndBulk <- merge(bulkPseudoCounts, sc_uc, by="row.names")
row.names(pseudoCountsScAndBulk) = pseudoCountsScAndBulk$Row.names
pseudoCountsScAndBulk = pseudoCountsScAndBulk[,-1]
pseudoCountsScAndBulk = tpmScAndBulk

tpmScAndBulk = MakeTPM(pseudoCountsScAndBulk)
colSums(tpmScAndBulk) #test


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
desm <- read.xlsx(paste0(folder, "/DesignMatrix.xlsx"), sheetName = "DesignMatrix")
cellTypes = as.numeric(desm[9, 2:106])
labs = as.numeric(desm[2, 2:106])
subCellTypes = as.numeric(desm[7, 2:106])
tissues = as.numeric(desm[3, 2:106])

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



###########################
## Figure 1
###########################

#for comparison:
#plotRLE(as.matrix(bulk_counts), outline=FALSE, ylim=c(-4, 4), col=as.numeric(cellTypes))
#plotRLE(as.matrix(bulk_tpm), outline=FALSE, ylim=c(-4, 4), col=as.numeric(cellTypes))
#plotRLE(as.matrix(bulk_tmm), outline=FALSE, ylim=c(-4, 4), col=as.numeric(cellTypes))

#so, instead of plotRLE we want to plot this using ggplot so we can put all of them together in a figure
#first: prepare data

library("matrixStats")

genRLEData <- function(df) {
  df = as.matrix(df)
  med = rowMedians(df)
  output = df
  for (i in 1:dim(output)[2]) {
    output[,i] = log2((df[,i] + 1) / (med + 1))
    
  }
  return (output)
}

a = genRLEData(tpmScAndBulk)
b = genRLEData(tmmScAndBulk)
c = genRLEData(quantileScAndBulk)

library(dplyr)
library(tidyr)

labvals = c(3,2,4,5,1,6,7,8,9,10)
labnames = c("Bulk 1", "Bulk 2", "Bulk 3", "Bulk 4", "Bulk 5", "SC HCA CB", "SC LC", "SC PBMC68k", "SC Mixed 10x data", "SC Melanoma")

plotSampleOrder = c(2:7,1,8:12,25,13:24,26:64,65:74,75:90,91:98,104:105,101:103,99:100)
plotSampleOrderFac = factor(1:length(plotSampleOrder), levels=plotSampleOrder)

#now with ggplot:
dfTPM <- a %>% 
  as.data.frame %>% 
  gather(., key= sample, value="value") %>%
  mutate(Dataset = rep(factor(labs, labvals, labnames),each=dim(a)[1]), 
         plotSampleOrder = rep(plotSampleOrderFac, each=dim(a)[1]))

dfTMM <- b %>% 
  as.data.frame %>% 
  gather(., key= sample, value="value") %>%
  mutate(Dataset = rep(factor(labs, labvals, labnames),each=dim(b)[1]), 
         plotSampleOrder = rep(plotSampleOrderFac, each=dim(b)[1]))

dfQ <- c %>% 
  as.data.frame %>% 
  gather(., key= sample, value="value") %>%
  mutate(Dataset = rep(factor(labs, labvals, labnames),each=dim(c)[1]), 
         plotSampleOrder = rep(plotSampleOrderFac, each=dim(c)[1]))


xline = c(0,dim(a)[2])
yline = c(0,0)
dfline = data.frame(xline,yline)

#Combine all three normalization methods:
df = rbind(dfTPM,dfTMM, dfQ)
len = dim(dfTPM)[1]
nm = rep(c(1,2,3), each=len)
nm = factor(nm, c(1,2,3), c("TPM/CPM", "TMM", "Quantile"))
df = cbind(df,nm)

ggplot(data= df , aes(x=plotSampleOrder, y=value, fill=Dataset)) +
  geom_boxplot(outlier.shape = NA, coef = 0) +
  labs(title="Relative Log Expression After Normalization", y="Relative Log Expression", x="Samples") +
  coord_cartesian(ylim=c(-2, 2)) +
  theme( axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
  geom_line(data = dfline, aes(x=xline, y=yline, fill=NA)) + 
  facet_wrap( ~ nm, nrow=3) 

###########################
## Figure 2 - PCA plots
###########################

#log transform the data and run PCA

a2 = prcomp(t(log2(tpmScAndBulk + 1)))
b2 = prcomp(t(log2(tmmScAndBulk + 1)))
c2 = prcomp(t(log2(quantileScAndBulk + 1)))
d2 = prcomp(t(log2(bcScAndBulk + 1)))

#calculate the explained variance for principle component 1 and 2 for each 
explVar <- function(pcaData, pc) {
  eigs <- pcaData$sdev^2
  return (eigs[pc] / sum(eigs))
}

x=c(1,1,1,1)
explPC1=c(explVar(a2, 1), explVar(b2, 1), explVar(c2, 1), explVar(d2, 1))
explPC2=c(explVar(a2, 2), explVar(b2, 2), explVar(c2, 2), explVar(d2, 2))
textPC1=c(paste0(round(100*explPC1[1],2), "% - Technical"),
          paste0(round(100*explPC1[2],2), "% - Technical"),
          paste0(round(100*explPC1[3],2), "% - Technical"),
          paste0(round(100*explPC1[4],2), "% - Cell type"))
textPC2=c(paste0(round(100*explPC2[1],2), "% - Cell type"),
          paste0(round(100*explPC2[2],2), "% - Cell type"),
          paste0(round(100*explPC2[3],2), "% - Cell type"),
          paste0(round(100*explPC2[4],2), "% - Unknown"))
nm = factor(c(1,2,3,4), c(1,2,3,4), c("TPM/CPM", "TMM", "Quantile", "TMM + ComBat"))
dat2 = data.frame(x,explPC1,explPC2,textPC1,nm)

summary(d2) #check that the PC1 explained variance fraction was correctly calculated -> OK!


#scale all PCs to be between 0 and 1
scalePCs <- function(d) {
  for (i in 1:(dim(d)[2])) {
    minval = min(d[,i])
    maxval = max(d[,i])
    d[,i] = (d[,i] - minval) / (maxval - minval)
  }
  return (d)
}
a2$x = scalePCs(a2$x)
b2$x = scalePCs(b2$x)
c2$x = scalePCs(c2$x)
d2$x = scalePCs(d2$x)
#flip PC2 upside down for c so the graphs are more similar. Same for PC1 for b
c2$x[,2] = -c2$x[,2] + 1
b2$x[,1] = -b2$x[,1] + 1




#plot(a2$x[,1], a2$x[,2])

dfTPM <- as.data.frame(a2$x)
dfTPM$labs <- factor(labs, labvals, labnames)
dfTPM$bulkAndCt <- factor((cellTypes - 1) + (labs > 5) * 2, c(0,1,2,3), c("Bulk B Cell", "Bulk T Cell", "Sc B Cell", "Sc T Cell"))
dfTPM$pc1Expl = explVar(a2, 1)
dfTPM$pc2Expl = explVar(a2, 2)

dfTMM <- as.data.frame(b2$x)
dfTMM$labs <- factor(labs, labvals, labnames)
dfTMM$bulkAndCt <- factor((cellTypes - 1) + (labs > 5) * 2, c(0,1,2,3), c("Bulk B Cell", "Bulk T Cell", "Sc B Cell", "Sc T Cell"))
dfTMM$pc1Expl = explVar(b2, 1)
dfTMM$pc2Expl = explVar(b2, 2)

dfQ <- as.data.frame(c2$x)
dfQ$labs <- factor(labs, labvals, labnames)
dfQ$bulkAndCt <- factor((cellTypes - 1) + (labs > 5) * 2, c(0,1,2,3), c("Bulk B Cell", "Bulk T Cell", "Sc B Cell", "Sc T Cell"))
dfQ$pc1Expl = explVar(c2, 1)
dfQ$pc2Expl = explVar(c2, 2)

dfBC <- as.data.frame(d2$x)
dfBC$labs <- factor(labs, labvals, labnames)
dfBC$bulkAndCt <- factor((cellTypes - 1) + (labs > 5) * 2, c(0,1,2,3), c("Bulk B Cell", "Bulk T Cell", "Sc B Cell", "Sc T Cell"))
dfBC$pc1Expl = explVar(d2, 1)
dfBC$pc2Expl = explVar(d2, 2)


#Combine all methods:
df = rbind(dfTPM,dfTMM, dfQ, dfBC)
len = dim(dfTPM)[1]
nm = rep(c(1,2,3,4), each=len)
nm = factor(nm, c(1,2,3,4), c("TPM/CPM", "TMM", "Quantile", "TMM + ComBat"))
df = cbind(df,nm)


p<-ggplot(df,aes(x=PC1,y=PC2,color=labs, shape=bulkAndCt))
p<-p + geom_point()
p<-p + labs(title="Visualization of Batch Effects")
p$labels$shape = "Sample properties"
p$labels$colour = "Dataset"
p<-p + coord_cartesian(xlim=c(-0.09, 1), ylim=c(-0.09, 1))#create room for PC label
p<-p + theme( axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) 

p<-p + geom_rect(data=dat2, mapping=aes(xmin=-0.114, xmax=-0.04, ymin=0, ymax=explPC2), color=NA, fill="blue", alpha=0.2, inherit.aes = FALSE)
p<-p + geom_rect(data=dat2, mapping=aes(xmin=-0.114, xmax=-0.04, ymin=explPC2, ymax=1), color=NA, alpha=0.0, inherit.aes = FALSE)
p<-p + geom_rect(data=dat2, mapping=aes(xmin=-0.114, xmax=-0.04, ymin=0, ymax=1), color="blue", fill=NA, inherit.aes = FALSE)

p<-p + geom_rect(data=dat2, mapping=aes(ymin=-0.12, ymax=-0.04, xmin=0, xmax=explPC1), color=NA, fill="blue", alpha=0.2, inherit.aes = FALSE)
p<-p + geom_rect(data=dat2, mapping=aes(ymin=-0.12, ymax=-0.04, xmin=explPC1, xmax=1), color=NA, alpha=0.0, inherit.aes = FALSE)
p<-p + geom_rect(data=dat2, mapping=aes(ymin=-0.12, ymax=-0.04, xmin=0, xmax=1), color="blue", fill=NA, inherit.aes = FALSE)

p<-p + geom_text(data=dat2, mapping=aes(y=-0.075, x=0.5, label=textPC1), size=3.8, color="black", inherit.aes = FALSE)
p<-p + geom_text(data=dat2, mapping=aes(x=-0.085, y=0.5, label=textPC2), size=3.8, color="black", angle = 90, inherit.aes = FALSE)
p<-p + facet_wrap( ~ nm, nrow=2) 
p

###########################
## Figure 3 first half - Importance of Factors
###########################

addon = 0.05
a3 = log2((tmmScAndBulk[,2] + addon)/(tmmScAndBulk[,3] + addon)) #technical replicates
b3 = log2((tmmScAndBulk[,13] + addon)/(tmmScAndBulk[,14] + addon)) #same lab, same individual
c3 = log2((tmmScAndBulk[,10] + addon)/(tmmScAndBulk[,11] + addon)) #same lab
d3 = log2((tmmScAndBulk[,2] + addon)/(tmmScAndBulk[,14] + addon)) #different lab
e3 = log2((tmmScAndBulk[,2] + addon)/(tmmScAndBulk[,8] + addon)) #different cell type 
f3 = log2((tmmScAndBulk[,8] + addon)/(tmmScAndBulk[,9] + addon)) #different sub cell type

lim = 4.5

a3 = a3[(a3 >= -lim) & (a3 <= lim)]
b3 = b3[(b3 >= -lim) & (b3 <= lim)]
c3 = c3[(c3 >= -lim) & (c3 <= lim)]
d3 = d3[(d3 >= -lim) & (d3 <= lim)]
e3 = e3[(e3 >= -lim) & (e3 <= lim)]
f3 = f3[(f3 >= -lim) & (f3 <= lim)]

breaks = seq(-4.55,4.55,by=0.1)

ha = hist(a3,breaks,plot=F)  #returns absolute  frequency
ax = ha$mids
ay = ha$counts
hb = hist(b3,breaks,plot=F)  #returns absolute  frequency
bx = hb$mids
by = hb$counts
hc = hist(c3,breaks,plot=F)  #returns absolute  frequency
cx = hc$mids
cy = hc$counts
hd = hist(d3,breaks,plot=F)  #returns absolute  frequency
dx = hd$mids
dy = hd$counts
he = hist(e3,breaks,plot=F)  #returns absolute  frequency
ex = he$mids
ey = he$counts
hf = hist(f3,breaks,plot=F)  #returns absolute  frequency
fx = hf$mids
fy = hf$counts

#x = c(ax,bx,cx,dx,ex,fx)
#y = c(ay,by,cy,dy,ey,fy)
#comparisons = rep(c(1,2,3,4,5,6),each=length(ax))
#comparisons = factor(comparisons, c(1,2,3,4,5,6), c("Techn. repl.", "Same pat.", "Diff. pat.", "Diff. lab", "Diff. cell type", "Diff sub cell type" ))

#skip technical replicates and same individual, too many curves in the figure
x = c(cx,dx,ex,fx)
y = c(cy,dy,ey,fy)
comparisons = rep(c(1,2,3,4),each=length(ax))
comparisons = factor(comparisons, c(1,2,3,4), c("Same lab and ct", "Diff. lab", "Diff. cell type", "Diff. sub cell type" ))


df = data.frame(x,y, comparisons)
#ggplot(data=df, aes(df$a)) +
#  geom_histogram(binwidth = 0.05) +
#  geom_histogram(aes = b, binwidth = 0.05)

#num = length(a3);
#sel = c(rep(0,num),rep(1,length(b)));


ggplot(data=df, aes(x=x, y=y, group=comparisons)) +
  geom_line(aes(colour=comparisons)) +
  labs(title="Histogram for Gene Log Fold Changes Between Samples", y="Number of genes", x="Log2 fold change")

#ggplot(data=df, aes(x=c)) +
#  geom_histogram(data=subset(df,sel == 0),fill = "red", alpha = 0.2, binwidth = 0.05) +
#  geom_histogram(data=subset(df,sel == 1),fill = "blue", alpha = 0.2, binwidth = 0.05) +
#  coord_cartesian(xlim=c(-5, 5))



#nonKallistoData = read.csv("C:/Work/MatlabCode/components/SCLib/ImportableData/tcellProfilesFromMatlab.txt", header=TRUE, row.names=1, sep = "\t", quote = "\"'", check.names=FALSE);
#colnames(nonKallistoData)
#totjoin = inner_join(rownames_to_column(as.data.frame(tmmScAndBulkNonFilt)), rownames_to_column(as.data.frame(nonKallistoData)))
#row.names(totjoin) = totjoin$rowname
#totjoin = totjoin[,-1]
#totjoin = MakeTPM(totjoin)
nonKallistoData = read.table(paste0(dataFolder, "AuthorsProcBulk.txt"), header = T,sep="\t")
#reorder for the plot
nonKallistoData = nonKallistoData[,c(3,1,2)]
kallistoData = tpmScAndBulkNonFilt[,c(19,33,12)]

compData = inner_join(rownames_to_column(as.data.frame(kallistoData)), rownames_to_column(as.data.frame(nonKallistoData)))
row.names(compData) = compData$rowname
compData = compData[,-1]
compData = MakeTPM(compData)



#colnames(totjoin[,c(19,156,33,141,12,167)])
#colnames(totjoin)

t = log2(compData[,1] + 0.05) #CD8%2b%20T%20Cells%20%28pluriselect%29%2c%20donor090309%2c%20donation1.CNhs12176.12186-129A8
t2 = log2(compData[,4] + 0.05) # CD8%2b%20T%20Cells%20%28pluriselect%29%2c%20donor090309%2c%20donation1.CNhs12176.12186-129A8
t3 = log2(compData[,2] + 0.05)# C001FRB1
t4 = log2(compData[,5] + 0.05)# C001FRB1
t5 = log2(compData[,3] + 0.05) #ENCFF088DIY
t6 = log2(compData[,6] + 0.05) #ENCFF088DIY


ds = rep(c(1,2,3,4,5,6),each=length(t))
Technology = factor(ds,c(1,2,3,4,5,6), c("Fantom5 Kallisto","Fantom5 Align", "BLUEPRINT Kallisto","BLUEPRINT Align", "ENCODE Kallisto", "ENCODE Align"))
xxxx = c(t,t2,t3,t4,t5,t6)
s = rep(rep(c(1,2), each=length(t)),3)
style = factor(s,c(1,2),c("Kallisto", "Align"))
Sample = factor(rep(c(1,2,3),each=length(t)*2), c(1,2,3), c("Fantom5 T Cells", "BLUEPRINT T Cells", "ENCODE T Cells"))
dftest = data.frame(xxxx,Technology, style, Sample)

xxxx = c(t,t2,t3,t4,t5,t6)
library("metagen")
lcols3 = cbbPalette[c(1,1,2,2,3,3)]

ggplot(data=dftest, aes(x=xxxx, group=Technology)) +
  geom_density(aes(colour=Technology, linetype = style)) +
  labs(title="Histogram for Gene Expression - Kallisto vs Alignment", y="Gene density", x="Log2(TPM + 0.05)") +
  coord_cartesian(xlim=c(-3.5, 15)) +#get rid of empty space to the left
  scale_color_manual(values=lcols3) +
  scale_linetype_manual(values=c("solid", "dashed",  "solid", "dashed", "solid", "dashed"))


plot(t,t2)
plot(t5,t6)

cor(t,t2)
cor(t3,t4)
cor(t5,t6)

sel = t >= 1
cor(t[sel],t2[sel])



hist(t6[t6 < 0],100)
plot(t,t2)
plot(t3,t4)
plot(t5,t6)

#Export genes to Broad:
#lowAlign = t2 < -2.5
#highK = t > 1
#sum(lowAlign & highK)
#hist(t[lowAlign & highK])
#strangeGenes = t[lowAlign & highK]
#names(strangeGenes) = row.names(totjoin[lowAlign & highK,])
#setwd("C:/Work/R/RNASeqCTProfileEval")
#write.table(strangeGenes, "genesHigherInKallisto.txt")




###########################
## Figure 4 - Lowly expressed genes are more lowly expressed in single-cell data - to be skipped!
###########################

#First show that there are more lowly expressed genes in 10x data
#test density plot
library("ggpubr")



#t = log2(tpmScAndBulkNonFilt[,1] + 0.05)
t = log2(tpmScAndBulkNonFilt[,8] + 0.05)
#t2 = log2(tpmScAndBulkNonFilt[,83] + 0.05)
#t2 = log2(rowMeans(tpmScAndBulkNonFilt[,83:90]) + 0.05)
t2 = log2(tpmScAndBulkNonFilt[,75] + 0.05)
#t2 = log2(rowMeans(tpmScAndBulkNonFilt[,c(75:82,101,103)]) + 0.05)
t3 = log2(tpmScAndBulkNonFilt[,104] + 0.05)#melanoma T cells

#testds = inner_join(rownames_to_column(as.data.frame(HCASCE_m10x)), rownames_to_column(as.data.frame(tpmScAndBulkNonFilt)))
#row.names(testds) = testds$rowname
#testds = testds[,-1]


#test with the test dataset
#t = log2(HCASCE_m10x[,1] + 0.05)
#t2 = log2(HCASCE_m10x[,3] + 0.05)
#t3 = log2(HCASCE_mc[,3] + 0.05)

#t = log2(testds[,1] + 0.05)
#t2 = log2(testds[,3] + 0.05)
#t3 = log2(testds[,5] + 0.05)


ds = rep(c(1,2,3),each=length(t))
Technology = factor(ds,c(1,2,3), c("Bulk","10x", "Smart-Seq2"))
xxxx = c(t,t2,t3)
dftest = data.frame(xxxx,Technology)
#ggdensity(dftest,x="xxxx", fill="Technology", color="dsf", rug=T )
#ggdensity(dftest,x="xxxx", fill=NA, color="Technology", rug=T ) + 
#  labs(title="Gene Expression Histograms per Technology", y="Gene density", x="Gene expression")

xxxx = c(t,t2,t3)

ggplot(data=dftest, aes(x=xxxx, group=Technology)) +
  geom_density(aes(colour=Technology)) +
  labs(title="Histogram for Gene Expression Across Technologies", y="Gene density", x="Log2(TPM + 0.05)") +
  coord_cartesian(xlim=c(-3.5, 15))#get rid of empty space to the left




#see what the lowly expressed genes in single-cell looks like in bulk:
#sel = EBVMerge[,1] < 0.25
#t = log2(EBVMerge[sel,2] + 0.05) # bulk
#t2 = log2(EBVMerge[sel,1] + 0.05) # single-cell
#ds = rep(c(1,2),each=length(t))
#dsf = factor(ds,c(1,2), c("Bulk","10x"))
#xxxx = c(t,t2)
#dftest = data.frame(xxxx,dsf)
#ggdensity(dftest,x="xxxx", fill="dsf", color="dsf", rug=T )
#hist(t,100)
#hist(t2,100)










#only look at the b cells
bcells = cellTypes == 1
labsb = labs[bcells]
ds = (labsb > 5) & (labsb != 10) #dropseq, 10x
ss = labsb == 10
bulk = (labsb <= 5)


bCellsTPM = tmmScAndBulkNonFilt[,bcells]
bBTPMMean = rowMeans(bCellsTPM)
dsBCellsTPM = bCellsTPM[,ds]
dsBTPMMean = rowMeans(dsBCellsTPM)
#have to filter very lowly expressed genes to get combat to work later
bCellsTPM = bCellsTPM[bBTPMMean > 0.1,]

bCellsTMM = TMMNorm(bCellsTPM)
bulkBCellsTMM = bCellsTMM[,bulk]
bBTMMMean = rowMeans(bulkBCellsTMM)
dsBCellsTMM = bCellsTMM[,ds]
dsBTMMMean = rowMeans(dsBCellsTMM)
ssBTMMMean = bCellsTMM[,ss] # there's really just one here...

dsTMMData = GetScBias(dsBTMMMean, bBTMMMean)
ssTMMData = GetScBias(ssBTMMMean, bBTMMMean)


#do quantile normalization on the bcells only
bCellsQ = normalize.quantiles(as.matrix(bCellsTPM))

bulkBCellsQ = bCellsQ[,bulk]
bBQMean = rowMeans(bulkBCellsQ)
dsBCellsQ = bCellsQ[,ds]
dsQMean = rowMeans(dsBCellsQ)
ssQMean = bCellsQ[,ss] # there's really just one here...

dsQData = GetScBias(dsQMean, bBQMean)
ssQData = GetScBias(ssQMean, bBQMean)

#do batch correction on the b cells only; skip lab 2,8-10 since we only have one sample each
d = as.matrix(log2(bCellsQ + 0.05))
e = as.matrix(log2(bCellsTMM + 0.05))
f = as.matrix(log2(bCellsTPM + 0.05))
sel = (labsb < 8) & (labsb != 2)
d = d[,sel]
e = e[,sel]
f = f[,sel]

batch = labsb[sel];

modcombat = model.matrix(~1, as.data.frame(batch))
combat_d = ComBat(dat=d, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
combat_e = ComBat(dat=e, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
combat_f = ComBat(dat=f, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

#transform back
bcTPMBCells = combat_f^2 - 0.05
bcTPMBCells[bcTPMBCells<0] = 0

bulkBCellsBCTPM = bcTPMBCells[,bulk[sel]]
bBBCTPMMean = rowMeans(bulkBCellsBCTPM)
dsBCellsBCTPM = bcTPMBCells[,ds[sel]]
dsBCTPMMean = rowMeans(dsBCellsBCTPM)

dsBCTPMData = GetScBias(dsBCTPMMean, bBBCTPMMean)


bcTMMBCells = combat_e^2 - 0.05
bcTMMBCells[bcTMMBCells<0] = 0

bulkBCellsBCTMM = bcTMMBCells[,bulk[sel]]
bBBCTMMMean = rowMeans(bulkBCellsBCTMM)
dsBCellsBCTMM = bcTMMBCells[,ds[sel]]
dsBCTMMMean = rowMeans(dsBCellsBCTMM)

dsBCTMMData = GetScBias(dsBCTMMMean, bBBCTMMMean)

bcQBCells = combat_edata^2 - 0.05
bcQBCells[bcQBCells<0] = 0

bulkBCellsBCQ = bcQBCells[,bulk[sel]]
bBBCQMean = rowMeans(bulkBCellsBCQ)
dsBCellsBCQ = bcQBCells[,ds[sel]]
dsBCQMean = rowMeans(dsBCellsBCQ)

dsBCQData = GetScBias(dsBCQMean, bBBCQMean)
dsBCQData = GetScBias(dsBCQMean * 10^6 / sum(dsBCQMean), bBBCQMean * 10^6 / sum(bBBCQMean))


#merge into a single df

df = as.data.frame(rbind(dsTMMData, ssTMMData, dsQData, ssQData, dsBCTPMData, dsBCTMMData, dsBCQData))
g = rep(c(1,2,3,4,5,6,7), each=dim(dsTMMData)[1])
Group = factor(g, c(1,2,3,4,5,6,7), c("TMM, 10x", "TMM, Mod. SmartSeq2", "Quantile, 10x", "Quantile, Mod. SmartSeq2    ", "TPM, 10x, ComBat", "TMM, 10x, ComBat", "Quantile, 10x, ComBat"))
color = Group
color[color == "TMM, 10x, ComBat"] = "TMM, 10x"
color[color == "Quantile, 10x, ComBat"] = "Quantile, 10x"
linetype = as.factor(c(rep(1, 4*dim(dsTMMData)[1]), rep(2, 3*dim(dsTMMData)[1])))

df = cbind(df,Group, color, linetype)
#remove all rows with NA
df = df[(!is.na(df[,"x"])) & (!is.na(df[,"y"])),]

#install.packages("metagen")
library("metagen")
lcols2 = cbbPalette[c(1,2,3,4,7,1,3)]

xline = c(-3,12)
yline = c(0,0)
dfline = data.frame(xline,yline)


plotDifferentDs = ggplot(data=df, aes(x=x, y=y, group=Group, color=Group)) +
  geom_line(data = dfline, inherit.aes=F, aes(x=xline, y=yline)) +
  geom_line(aes(linetype = Group)) +
  labs(title="Single-Cell Biases vs Bulk - Comparison Across Datasets", y="Relative Log Expression", x="Gene Expression (pseudo-TPM)") +
  scale_color_manual(values=lcols2) +
  scale_linetype_manual(values=c("solid", "solid",  "solid", "solid", "dashed", "dashed", "dashed"))

#now create a plot from the HCA single-cell evaluation dataset
dsccTPMData = GetScBias(HCASCE_m10x[,5], HCASCE_m10x[,1], -3, 11)
dscuTPMData = GetScBias(HCASCE_m10x[,3], HCASCE_m10x[,1], -3, 11)
dsccTMMData = GetScBias(HCASCE_m10xn[,5], HCASCE_m10xn[,1], -3, 11)
dscuTMMData = GetScBias(HCASCE_m10xn[,3], HCASCE_m10xn[,1], -3, 11)
ssccTPMData = GetScBias(HCASCE_mc[,3], HCASCE_mc[,1], -3, 11)
ssccTMMData = GetScBias(HCASCE_mcn[,3], HCASCE_mcn[,1], -3, 11)

dfeval = as.data.frame(rbind(dsccTPMData, dscuTPMData, dsccTMMData, dscuTMMData, ssccTPMData, ssccTMMData))
g = rep(c(1,2,3,4,5,6), each=dim(dscuTPMData)[1])
Sample = factor(g, c(1,2,3,4,5,6), c("TPM, 10x, counts", "TPM, 10x, UMIs", "TMM, 10x, counts", "TMM, 10x, UMIs", "TPM, Smart-Seq2, counts", "TMM, Smart-Seq2, counts"))
#color = Sample
#color[color == "TPM, Smart-Seq2, cortex 1"] = "TPM, 10x, cortex 1"
#color[color == "TPM, Smart-Seq2, cortex 2"] = "TPM, 10x, cortex 2"
#color[color == "TMM, Smart-Seq2, cortex 1"] = "TMM, 10x, cortex 1"
#color[color == "TMM, Smart-Seq2, cortex 2"] = "TMM, Smart-Seq2, cortex 2"
#linetype = as.factor(rep(c(1,2), each = dim(dsTMMData)[1], times=3))

dfeval = cbind(dfeval,Sample)
#remove all rows with NA
dfeval = dfeval[(!is.na(dfeval[,"x"])) & (!is.na(dfeval[,"y"])),]

#install.packages("metagen")
library("metagen")
lcols = cbbPalette[c(1,1,2,2,3,4)]

xline = c(-3,12)
yline = c(0,0)
dfline = data.frame(xline,yline)


plotEval = ggplot(data=dfeval, aes(x=x, y=y, group=Sample, color=Sample)) +
  geom_line(data = dfline, inherit.aes=F, aes(x=xline, y=yline)) +
  geom_line(aes(linetype = Sample)) +
  labs(title="Single-Cell Biases vs Bulk - Different Technologies on the Same Samples", y="Relative Log Expression", x="Bulk Gene Expression (pseudo-TPM)") +
  scale_color_manual(values=lcols) +
  scale_linetype_manual(values=c("solid", "dashed",  "solid", "dashed", "solid", "solid"))


ggarrange(
  plotEval,                # First row with line plot
  # Second row with box and dot plots
  plotDifferentDs, 
  nrow = 2, 
  labels = c("A","B")       # Label of the line plot
) 


t = log2(HCASCE_m10x[,3] + 0.05) #sc UMIs
t2 = log2(HCASCE_m10x[,1] + 0.05) #bulk
t3 = log2(HCASCE_m10x[,5] + 0.05) #sc counts


ds = rep(c(1,2,3),each=length(t))
Technology = factor(ds,c(1,2,3), c("sc UMIs","bulk","sc counts"))
xxxx = c(t,t2,t3)
s = rep(rep(c(1,2,3), each=length(t)))
dftest = data.frame(xxxx,Technology)

ggplot(data=dftest, aes(x=xxxx, group=Technology)) +
  geom_density(aes(colour=Technology)) +
  labs(title="Histogram for Gene Expression - UMI vs bulk", y="Gene density", x="Log2(TPM + 0.05)") +
  coord_cartesian(xlim=c(-3.5, 15)) #get rid of empty space to the left

sel = t3 < 0
hist(t3[t3 < 0],50)


###########################
## Figure 5 - Technical biases between single-cell and bulk
###########################

#Four things: 1. UMIs vs counts (removed counts' fraction), 2. gene length, 3. GC content 4. GC content tail

cort1 = read.table(paste0(dataFolder, "/ScVsBulkCortex1.txt"), sep="\t")
cort2 = read.table(paste0(dataFolder, "/ScVsBulkCortex2.txt"), sep="\t")

#check correlation between UMIFrac for cort1 and cort2. Not super, but it seems it is there
plot(cort1$remUMIFrac, cort2$remUMIFrac)

#Remove all lowly expressed genes, so much noise there
#filt = (cort1$logUMITMM >= log2(1.05)) & (cort1$logBulkTMM >= log2(1.05))   #filter on pseudo-TPM of 1
#filtCortExpr1 = cort1[filt,]
#filt2 = (cort2$logUMITMM >= log2(1.05)) & (cort2$logBulkTMM >= log2(1.05)) #filter on pseudo-TPM of 1
#filtCortExpr2 = cort2[filt2,]

#Don't filter on gene expression, it is very likely that this will introduce some bias!
filtCortExpr1 = cort1
filtCortExpr2 = cort2



#Remove outliers (quite few genes)
filtGeneLength1 = filtCortExpr1$geneLength < 10000
filtGeneLength2 = filtCortExpr2$geneLength < 10000
filtGCFull1 = (filtCortExpr1$gcFullLength >= 0.33) & (filtCortExpr1$gcFullLength <= 0.67)
filtGCFull2 = (filtCortExpr2$gcFullLength >= 0.33) & (filtCortExpr2$gcFullLength <= 0.67)
filtGCTail1 = (filtCortExpr1$gcTail >= 0.2) & (filtCortExpr1$gcTail <= 0.65)
filtGCTail2 = (filtCortExpr2$gcTail >= 0.2) & (filtCortExpr2$gcTail <= 0.65)
filtAll1 = filtGeneLength1 & filtGCFull1 & filtGCTail1
filtAll2 = filtGeneLength2 & filtGCFull2 & filtGCTail2
sum(!filtAll1) #number of discarded outliers cortex 1, 98
sum(!filtAll2) #number of discarded outliers cortex 2, 103
sum(filtAll1) #number of genes left cortex 1, 16510
sum(filtAll2) #number of genes left cortex 2, 16692

filtCort1 = filtCortExpr1[filtAll1,]
filtCort2 = filtCortExpr2[filtAll2,]


#1 Removed counts' fraction
###########################
resCort1RemUMIFrac = genScToBulkCovGraphs(filtCort1, logUMITMM ~ remUMIFracOtherSample, LogUMIDivBulk ~ remUMIFracOtherSample, 7, "UMI copy fraction")
resCort2RemUMIFrac = genScToBulkCovGraphs(filtCort2, logUMITMM ~ remUMIFracOtherSample, LogUMIDivBulk ~ remUMIFracOtherSample, 7, "UMI copy fraction")

#2. Gene length
###########################

#filter out outlier genes above 10000 in length
resCort1GeneLength = genScToBulkCovGraphs(filtCort1, logUMITMM ~ geneLength, LogUMIDivBulk ~ geneLength, 5, "Gene length")
resCort2GeneLength = genScToBulkCovGraphs(filtCort2, logUMITMM ~ geneLength, LogUMIDivBulk ~ geneLength, 5, "Gene length")

#3. GC Content full length
###########################
resCort1GCFullLength = genScToBulkCovGraphs(filtCort1, logUMITMM ~ gcFullLength, LogUMIDivBulk ~ gcFullLength, 3, "GC content entire transcript")
resCort2GCFullLength = genScToBulkCovGraphs(filtCort2, logUMITMM ~ gcFullLength, LogUMIDivBulk ~ gcFullLength, 3, "GC content entire transcript")

#4. GC Content tail
###########################
resCort1GCTail = genScToBulkCovGraphs(filtCort1, logUMITMM ~ gcTail, LogUMIDivBulk ~ gcTail, 4, "GC content transcript tail")
resCort2GCTail = genScToBulkCovGraphs(filtCort2, logUMITMM ~ gcTail, LogUMIDivBulk ~ gcTail, 4, "GC content transcript tail")

#5. Test to regress out all covariates:
###########################


#Regress out all covariates
formAll = LogUMIDivBulk ~ remUMIFracOtherSample + geneLength + gcFullLength + gcTail
loess_fitAll1 <- loess(formAll, filtCort1)
lmAll1 <- lm(formAll, filtCort1)
loess_fitAll2 <- loess(formAll, filtCort2)
lmAll2 <- lm(formAll, filtCort2)

cort1RegrAllLoess = regrOutUMIVsBulk(filtCort1, loess_fitAll1)
cort1RegrAllLin = regrOutUMIVsBulk(filtCort1, lmAll1)
cort2RegrAllLoess = regrOutUMIVsBulk(filtCort2, loess_fitAll2)
cort2RegrAllLin = regrOutUMIVsBulk(filtCort2, lmAll2)

corAllLoess1 = corrUMIVsBulk(cort1RegrAllLoess)
corAllLin1 = corrUMIVsBulk(cort1RegrAllLin)
corAllLoess2 = corrUMIVsBulk(cort2RegrAllLoess)
corAllLin2 = corrUMIVsBulk(cort2RegrAllLin)

#regress out all but gcTail
formAllButGCTail = LogUMIDivBulk ~ remUMIFracOtherSample + geneLength + gcFullLength
loess_fitAllButGCTail1 <- loess(formAllButGCTail, filtCort1)
lmAllButGCTail1 <- lm(formAllButGCTail, filtCort1)
loess_fitAllButGCTail2 <- loess(formAllButGCTail, filtCort2)
lmAllButGCTail2 <- lm(formAllButGCTail, filtCort2)

cort1RegrAllButGCTailLoess = regrOutUMIVsBulk(filtCort1, loess_fitAllButGCTail1)
cort1RegrAllButGCTailLin = regrOutUMIVsBulk(filtCort1, lmAllButGCTail1)
cort2RegrAllButGCTailLoess = regrOutUMIVsBulk(filtCort2, loess_fitAllButGCTail2)
cort2RegrAllButGCTailLin = regrOutUMIVsBulk(filtCort2, lmAllButGCTail2)

corAllButGCTailLoess1 = corrUMIVsBulk(cort1RegrAllButGCTailLoess)
corAllButGCTailLin1 = corrUMIVsBulk(cort1RegrAllButGCTailLin)
corAllButGCTailLoess2 = corrUMIVsBulk(cort2RegrAllButGCTailLoess)
corAllButGCTailLin2 = corrUMIVsBulk(cort2RegrAllButGCTailLin)

#regress out remUMIFrac and gcFullLength (i.e. remove gene length as well)
formRemUMIFracAndGCFullLength = LogUMIDivBulk ~ remUMIFracOtherSample + gcFullLength
loess_fitRemUMIFracAndGCFullLength1 <- loess(formRemUMIFracAndGCFullLength, filtCort1)
lmRemUMIFracAndGCFullLength1 <- lm(formRemUMIFracAndGCFullLength, filtCort1)
loess_fitRemUMIFracAndGCFullLength2 <- loess(formRemUMIFracAndGCFullLength, filtCort2)
lmRemUMIFracAndGCFullLength2 <- lm(formRemUMIFracAndGCFullLength, filtCort2)

cort1RegrRemUMIFracAndGCFullLengthLoess = regrOutUMIVsBulk(filtCort1, loess_fitRemUMIFracAndGCFullLength1)
cort1RegrRemUMIFracAndGCFullLengthLin = regrOutUMIVsBulk(filtCort1, lmRemUMIFracAndGCFullLength1)
cort2RegrRemUMIFracAndGCFullLengthLoess = regrOutUMIVsBulk(filtCort2, loess_fitRemUMIFracAndGCFullLength2)
cort2RegrRemUMIFracAndGCFullLengthLin = regrOutUMIVsBulk(filtCort2, lmRemUMIFracAndGCFullLength2)

corRemUMIFracAndGCFullLengthLoess1 = corrUMIVsBulk(cort1RegrRemUMIFracAndGCFullLengthLoess)
corRemUMIFracAndGCFullLengthLin1 = corrUMIVsBulk(cort1RegrRemUMIFracAndGCFullLengthLin)
corRemUMIFracAndGCFullLengthLoess2 = corrUMIVsBulk(cort2RegrRemUMIFracAndGCFullLengthLoess)
corRemUMIFracAndGCFullLengthLin2 = corrUMIVsBulk(cort2RegrRemUMIFracAndGCFullLengthLin)

formGeneLengthAndGCFullLength = LogUMIDivBulk ~ geneLength + gcFullLength
loess_fitGeneLengthAndGCFullLength1 <- loess(formGeneLengthAndGCFullLength, filtCort1)
cort1GeneLengthAndGCFullLengthLoess = regrOutUMIVsBulk(filtCort1, loess_fitGeneLengthAndGCFullLength1)
corGeneLengthAndGCFullLengthLoess1 = corrUMIVsBulk(cort1GeneLengthAndGCFullLengthLoess)


#now collect all the data
#create all columns
#indices in res variables, regressing out var: 3 = none, 4 = loess, 5 = linear
#rows: none, UMICF, Gene Length, GC Content, 
#      GC Content Tail, All, UMICF + GC content + Gene length, UMICF + GC content
c1Lin = c(resCort1RemUMIFrac[[3]], resCort1RemUMIFrac[[5]], resCort1GeneLength[[5]], resCort1GCFullLength[[5]],
          resCort1GCTail[[5]], corAllLin1, corAllButGCTailLin1, corRemUMIFracAndGCFullLengthLin1)
c2Lin = c(resCort2RemUMIFrac[[3]], resCort2RemUMIFrac[[5]], resCort2GeneLength[[5]], resCort2GCFullLength[[5]],
          resCort2GCTail[[5]], corAllLin2, corAllButGCTailLin2, corRemUMIFracAndGCFullLengthLin2)
meanLin = (c1Lin + c2Lin)/2
c1Loess = c(resCort1RemUMIFrac[[3]], resCort1RemUMIFrac[[4]], resCort1GeneLength[[4]], resCort1GCFullLength[[4]],
            resCort1GCTail[[4]], corAllLoess1, corAllButGCTailLoess1, corRemUMIFracAndGCFullLengthLoess1)
c2Loess = c(resCort2RemUMIFrac[[3]], resCort2RemUMIFrac[[4]], resCort2GeneLength[[4]], resCort2GCFullLength[[4]],
            resCort2GCTail[[4]], corAllLoess2, corAllButGCTailLoess2, corRemUMIFracAndGCFullLengthLoess2)
meanLoess = (c1Loess + c2Loess)/2

corTable = data.frame(c1Lin, c2Lin, meanLin, c1Loess, c2Loess, meanLoess)
round(corTable,3)

x = factor(1:8, 1:8, c("None", "UMICF", "Gene length", "GC content", "GC cont. tail", "UMICF\nGC content\ngene length\nGC cont. tail", "UMICF\nGC content\ngene length", "UMICF\nGC content"))
Fit = factor(rep(c(1,2), each=length(meanLin)), 1:2, c("linear", "loess"))
y = c(meanLin, meanLoess)

dfPlot = data.frame(x, y, Fit)

ggplot(data=dfPlot, aes(x=x, y=y, fill=Fit)) +
  geom_bar(stat="identity",position=position_dodge()) +
  coord_cartesian(ylim=c(0.8, 0.87)) +
  labs(title="Mean Effect of Regressing out Technical Covariates", y="Correlation, 10x vs bulk, log transformed data", x="Covariates regressed out") #+
  #theme(axis.text.x = element_text(angle=65, vjust=0.6))
  




  

#create a combined plot of all covariates:
###########################################
plots = list(resCort1RemUMIFrac[[2]], resCort1GeneLength[[2]], resCort1GCFullLength[[2]], resCort1GCTail[[2]])
library("ggpubr")

fig = ggarrange( #when exporting this, make the y size larger(y=620)
  plotlist = plots, 
  nrow = 2, 
  ncol = 2,
  labels = c("A","B","C","D") 
) 
 
annotate_figure(fig,
                top = text_grob("Log2 Fold Change Between 10x and Bulk vs Technical Covariates", face = "bold", size = 14))
fig


ind = sort(filtCort1$geneLength, index.return=T)
filtCort1Sort2 = filtCort1[ind$ix,];
#remove a few outlier genes, the graphs are zoomed out otherwise...
outliers1 = (filtCort1Sort2$geneLength > 10000)

filtCort1Sort2 = filtCort1Sort2[!outliers1,]

plot(filtCort1Sort2$geneLength, filtCort1Sort2$logUMITMM)
loess_fit <- loess(logUMITMM ~ geneLength, filtCort1Sort2)
lines(filtCort1Sort2$geneLength, predict(loess_fit), col = "red")

plot(filtCort1Sort2$geneLength, filtCort1Sort2$LogUMIDivBulk)
loess_fit4 <- loess(LogUMIDivBulk ~ geneLength, filtCort1Sort2)
lines(filtCort1Sort2$geneLength, predict(loess_fit4), col = "red")
#make linear fit as well
lm4 <- lm(LogUMIDivBulk ~ geneLength, filtCort1Sort2)
lines(filtCort1Sort2$geneLength, predict(lm4), col = "blue")

#Regress out gene length
cort1RegrGeneLength = regrOutUMIVsBulk(filtCort1Sort2, loess_fit4)
cort1RegrGeneLengthLin = regrOutUMIVsBulk(filtCort1Sort2, lm4)

corrUMIVsBulk(filtCort1Sort2)
corrUMIVsBulk(cort1RegrGeneLength)
corrUMIVsBulk(cort1RegrGeneLengthLin)
plot(cort1RegrUMI$geneLength, cort1RegrUMI$LogUMIDivBulk)

#3. GC Content full length
###########################

ind = sort(filtCort1$gcFullLength, index.return=T)
filtCort1Sort3 = filtCort1[ind$ix,];
#remove a few outlier genes, the graphs are zoomed out otherwise...
outliers2 = (filtCort1Sort3$gcFullLength < 0.33) | (filtCort1Sort3$gcFullLength > 0.67)

filtCort1Sort3 = filtCort1Sort3[!outliers2,]

plot(filtCort1Sort3$gcFullLength, filtCort1Sort3$logUMITMM)
loess_fit <- loess(logUMITMM ~ gcFullLength, filtCort1Sort3)
lines(filtCort1Sort3$gcFullLength, predict(loess_fit), col = "red")

plot(filtCort1Sort3$gcFullLength, filtCort1Sort3$LogUMIDivBulk)
loess_fit5 <- loess(LogUMIDivBulk ~ gcFullLength, filtCort1Sort3)
lines(filtCort1Sort3$gcFullLength, predict(loess_fit5), col = "red")
#make linear fit as well
lm5 <- lm(LogUMIDivBulk ~ gcFullLength, filtCort1Sort3)
lines(filtCort1Sort3$gcFullLength, predict(lm5), col = "blue")

#Regress out gc content
cort1RegrGcFullLength = regrOutUMIVsBulk(filtCort1Sort3, loess_fit5)
cort1RegrGcFullLengthLin = regrOutUMIVsBulk(filtCort1Sort3, lm5)

corrUMIVsBulk(filtCort1Sort3)
corrUMIVsBulk(cort1RegrGcFullLength)
corrUMIVsBulk(cort1RegrGcFullLengthLin)
plot(cort1RegrGcFullLength$gcFullLength, cort1RegrGcFullLength$LogUMIDivBulk)

#look at correlation between gc content and UMI rem fraction
plot(filtCort1$remUMIFrac, filtCort1$gcFullLength)
cor(filtCort1$remUMIFrac, filtCort1$gcFullLength)

#4. GC Content tail
###########################

ind = sort(filtCort1$gcTail, index.return=T)
filtCort1Sort4 = filtCort1[ind$ix,];
#remove a few outlier genes, the graphs are zoomed out otherwise...
outliers3 = (filtCort1Sort4$gcTail < 0.2) | (filtCort1Sort4$gcTail > 0.65)

filtCort1Sort4 = filtCort1Sort4[!outliers3,]

plot(filtCort1Sort4$gcTail, filtCort1Sort4$logUMITMM)
loess_fit <- loess(logUMITMM ~ gcTail, filtCort1Sort4)
lines(filtCort1Sort4$gcTail, predict(loess_fit), col = "red")

plot(filtCort1Sort4$gcTail, filtCort1Sort4$LogUMIDivBulk)
loess_fit6 <- loess(LogUMIDivBulk ~ gcTail, filtCort1Sort4)
lines(filtCort1Sort4$gcTail, predict(loess_fit6), col = "red")
#make linear fit as well
lm6 <- lm(LogUMIDivBulk ~ gcTail, filtCort1Sort4)
lines(filtCort1Sort4$gcTail, predict(lm6), col = "blue")

#Regress out gc content
cort1RegrGcTail = regrOutUMIVsBulk(filtCort1Sort4, loess_fit6)
cort1RegrGcTailLin = regrOutUMIVsBulk(filtCort1Sort4, lm6)

corrUMIVsBulk(filtCort1Sort4)
corrUMIVsBulk(cort1RegrGcTail)
corrUMIVsBulk(cort1RegrGcTailLin)
plot(cort1RegrGcTail$gcTail, cort1RegrGcTail$LogUMIDivBulk)

#look at correlation between gc content and UMI rem fraction
plot(filtCort1$remUMIFrac, filtCort1$gcTail)
cor(filtCort1$remUMIFrac, filtCort1$gcTail)

#look at correlation between tail gc content and full gc content
plot(filtCort1$gcFullLength, filtCort1$gcTail)
cor(filtCort1$gcFullLength, filtCort1$gcTail)



#make a histogram plot displaying the improvement
c = c(cort1RegrUMI$LogUMIDivBulk, filtCort1Sort1$LogUMIDivBulk)
df9 = as.data.frame(c)
num = length(cort1RegrUMI$LogUMIDivBulk);
sel = c(rep(0,num),rep(1,num));
ggplot(data=df9, aes(x=c)) +
  geom_histogram(data=subset(df9,sel == 0),fill = "red", alpha = 0.2, binwidth = 0.2) +
  geom_histogram(data=subset(df9,sel == 1),fill = "blue", alpha = 0.2, binwidth = 0.2)



