library("ggplot2")
library("Seurat")

dataFolder = "C:/Work/R/RNASeqCTProfileEval/"
source(paste0(dataFolder, "FigureHelpFunc.R"))

###########################
## Figure 4 and 5 - Technical biases between single-cell and bulk
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
sum(!filtAll1) #number of discarded outliers cortex 1, 232
sum(!filtAll2) #number of discarded outliers cortex 2, 232
sum(filtAll1) #number of genes left cortex 1
sum(filtAll2) #number of genes left cortex 2


filtCort1 = filtCortExpr1[filtAll1,]
filtCort2 = filtCortExpr2[filtAll2,]

sum(is.na(filtCort1$remUMIFracOtherSample))
sum(is.na(filtCort2$remUMIFracOtherSample))


#experiment with not using umi copy fraction for genes with very few counts - it seems 5 UMIs or less gives only noise
filtCortRemRUFExpr1 = filtCortExpr1;
filtCortRemRUFExpr1$remUMIFracOtherSample[filtCortExpr2$UMI < 6] = NA
#filtCortRemRUFExpr1$remUMIFracOtherSample[filtCortExpr1$UMITPM < 0.001 | filtCortExpr1$bulkTPM < 0.001] = NA
filtCortRemRUF1 = filtCortRemRUFExpr1[filtAll1,]
#resCort1RemUMIFracSpec = genScToBulkCovGraphs(filtCortRemRUF1, logUMITMM ~ remUMIFracOtherSample, LogUMIDivBulk ~ remUMIFracOtherSample, 7, "UMI copy fraction")

filtCortRemRUFExpr2 = filtCortExpr2;
filtCortRemRUFExpr2$remUMIFracOtherSample[filtCortExpr1$UMI < 6] = NA
#filtCortRemRUFExpr1$remUMIFracOtherSample[filtCortExpr1$UMITPM < 0.001 | filtCortExpr1$bulkTPM < 0.001] = NA
filtCortRemRUF2 = filtCortRemRUFExpr2[filtAll2,]
#resCort2RemUMIFracSpec = genScToBulkCovGraphs(filtCortRemRUF2, logUMITMM ~ remUMIFracOtherSample, LogUMIDivBulk ~ remUMIFracOtherSample, 7, "UMI copy fraction")


#0 Copies per UMI
###########################
#resCort1CopiesPerUMI = genScToBulkCovGraphs(filtCort1, logUMITMM ~ CopiesPerUMIOtherSample, LogUMIDivBulk ~ CopiesPerUMIOtherSample, 14, "Copies per UMI")
#resCort2RemUMIFrac = genScToBulkCovGraphs(filtCort2, logUMITMM ~ remUMIFracOtherSample, LogUMIDivBulk ~ remUMIFracOtherSample, 7, "UMI copy fraction")


#1 Removed counts' fraction
###########################
#resCort1RemUMIFrac = genScToBulkCovGraphs(filtCort1, logUMITMM ~ remUMIFracOtherSample, LogUMIDivBulk ~ remUMIFracOtherSample, 7, "UMI copy fraction")
#resCort2RemUMIFrac = genScToBulkCovGraphs(filtCort2, logUMITMM ~ remUMIFracOtherSample, LogUMIDivBulk ~ remUMIFracOtherSample, 7, "UMI copy fraction")
#use the data where UMI copy fraction is removed where it is based on too few molecules
resCort1RemUMIFrac = genScToBulkCovGraphs(filtCortRemRUF1, logUMITMM ~ remUMIFracOtherSample, LogUMIDivBulk ~ remUMIFracOtherSample, 7, "UMI copy fraction")
resCort2RemUMIFrac = genScToBulkCovGraphs(filtCortRemRUF2, logUMITMM ~ remUMIFracOtherSample, LogUMIDivBulk ~ remUMIFracOtherSample, 7, "UMI copy fraction")

#2. Gene length
###########################

#filter out outlier genes above 10000 in length
resCort1GeneLength = genScToBulkCovGraphs(filtCortRemRUF1, logUMITMM ~ geneLength, LogUMIDivBulk ~ geneLength, 5, "Gene length")
resCort2GeneLength = genScToBulkCovGraphs(filtCortRemRUF2, logUMITMM ~ geneLength, LogUMIDivBulk ~ geneLength, 5, "Gene length")

#3. GC Content full length
###########################
resCort1GCFullLength = genScToBulkCovGraphs(filtCortRemRUF1, logUMITMM ~ gcFullLength, LogUMIDivBulk ~ gcFullLength, 3, "GC content entire transcript")
resCort2GCFullLength = genScToBulkCovGraphs(filtCortRemRUF2, logUMITMM ~ gcFullLength, LogUMIDivBulk ~ gcFullLength, 3, "GC content entire transcript")

#4. GC Content tail
###########################
resCort1GCTail = genScToBulkCovGraphs(filtCortRemRUF1, logUMITMM ~ gcTail, LogUMIDivBulk ~ gcTail, 4, "GC content transcript tail")
resCort2GCTail = genScToBulkCovGraphs(filtCortRemRUF2, logUMITMM ~ gcTail, LogUMIDivBulk ~ gcTail, 4, "GC content transcript tail")

#5. Test to regress out all covariates:
###########################


#Regress out all covariates
formAll = LogUMIDivBulk ~ remUMIFracOtherSample + geneLength + gcFullLength + gcTail
loess_fitAll1 <- loess(formAll, filtCortRemRUF1)
lmAll1 <- lm(formAll, filtCortRemRUF1)
loess_fitAll2 <- loess(formAll, filtCortRemRUF2)
lmAll2 <- lm(formAll, filtCortRemRUF2)

cort1RegrAllLoess = regrOutUMIVsBulk(filtCortRemRUF1, loess_fitAll1)
cort1RegrAllLin = regrOutUMIVsBulk(filtCortRemRUF1, lmAll1)
cort2RegrAllLoess = regrOutUMIVsBulk(filtCortRemRUF2, loess_fitAll2)
cort2RegrAllLin = regrOutUMIVsBulk(filtCortRemRUF2, lmAll2)

corAllLoess1 = corrUMIVsBulk(cort1RegrAllLoess)
corAllLin1 = corrUMIVsBulk(cort1RegrAllLin)
corAllLoess2 = corrUMIVsBulk(cort2RegrAllLoess)
corAllLin2 = corrUMIVsBulk(cort2RegrAllLin)

#regress out all but gcTail
formAllButGCTail = LogUMIDivBulk ~ remUMIFracOtherSample + geneLength + gcFullLength
loess_fitAllButGCTail1 <- loess(formAllButGCTail, filtCortRemRUF1)
lmAllButGCTail1 <- lm(formAllButGCTail, filtCortRemRUF1)
loess_fitAllButGCTail2 <- loess(formAllButGCTail, filtCortRemRUF2)
lmAllButGCTail2 <- lm(formAllButGCTail, filtCortRemRUF2)

cort1RegrAllButGCTailLoess = regrOutUMIVsBulk(filtCortRemRUF1, loess_fitAllButGCTail1)
cort1RegrAllButGCTailLin = regrOutUMIVsBulk(filtCortRemRUF1, lmAllButGCTail1)
cort2RegrAllButGCTailLoess = regrOutUMIVsBulk(filtCortRemRUF2, loess_fitAllButGCTail2)
cort2RegrAllButGCTailLin = regrOutUMIVsBulk(filtCortRemRUF2, lmAllButGCTail2)

corAllButGCTailLoess1 = corrUMIVsBulk(cort1RegrAllButGCTailLoess)
corAllButGCTailLin1 = corrUMIVsBulk(cort1RegrAllButGCTailLin)
corAllButGCTailLoess2 = corrUMIVsBulk(cort2RegrAllButGCTailLoess)
corAllButGCTailLin2 = corrUMIVsBulk(cort2RegrAllButGCTailLin)

#regress out remUMIFrac and gcFullLength (i.e. remove gene length as well)
formRemUMIFracAndGCFullLength = LogUMIDivBulk ~ remUMIFracOtherSample + gcFullLength
loess_fitRemUMIFracAndGCFullLength1 <- loess(formRemUMIFracAndGCFullLength, filtCortRemRUF1)
lmRemUMIFracAndGCFullLength1 <- lm(formRemUMIFracAndGCFullLength, filtCortRemRUF1)
loess_fitRemUMIFracAndGCFullLength2 <- loess(formRemUMIFracAndGCFullLength, filtCortRemRUF2)
lmRemUMIFracAndGCFullLength2 <- lm(formRemUMIFracAndGCFullLength, filtCortRemRUF2)

cort1RegrRemUMIFracAndGCFullLengthLoess = regrOutUMIVsBulk(filtCortRemRUF1, loess_fitRemUMIFracAndGCFullLength1)
cort1RegrRemUMIFracAndGCFullLengthLin = regrOutUMIVsBulk(filtCortRemRUF1, lmRemUMIFracAndGCFullLength1)
cort2RegrRemUMIFracAndGCFullLengthLoess = regrOutUMIVsBulk(filtCortRemRUF2, loess_fitRemUMIFracAndGCFullLength2)
cort2RegrRemUMIFracAndGCFullLengthLin = regrOutUMIVsBulk(filtCortRemRUF2, lmRemUMIFracAndGCFullLength2)

corRemUMIFracAndGCFullLengthLoess1 = corrUMIVsBulk(cort1RegrRemUMIFracAndGCFullLengthLoess)
corRemUMIFracAndGCFullLengthLin1 = corrUMIVsBulk(cort1RegrRemUMIFracAndGCFullLengthLin)
corRemUMIFracAndGCFullLengthLoess2 = corrUMIVsBulk(cort2RegrRemUMIFracAndGCFullLengthLoess)
corRemUMIFracAndGCFullLengthLin2 = corrUMIVsBulk(cort2RegrRemUMIFracAndGCFullLengthLin)

formGeneLengthAndGCFullLength = LogUMIDivBulk ~ geneLength + gcFullLength
loess_fitGeneLengthAndGCFullLength1 <- loess(formGeneLengthAndGCFullLength, filtCortRemRUF1)
cort1GeneLengthAndGCFullLengthLoess = regrOutUMIVsBulk(filtCortRemRUF1, loess_fitGeneLengthAndGCFullLength1)
corGeneLengthAndGCFullLengthLoess1 = corrUMIVsBulk(cort1GeneLengthAndGCFullLengthLoess)

plot(cort1GeneLengthAndGCFullLengthLoess$logBulkTMM, cort1GeneLengthAndGCFullLengthLoess$logUMITMM)
plot(filtCortRemRUF1$logBulkTMM, filtCortRemRUF1$logUMITMM)

#check correlation between cortex 1 and 2 as comparison
print("Correlation values in the text");
ds2 = cbind(filtCortRemRUF1$logBulkTMM, filtCortRemRUF1$logUMITMM, filtCortRemRUF2$logBulkTMM, filtCortRemRUF2$logUMITMM)
cor(ds2[,1], ds2[,3])
cor(ds2[,2], ds2[,4])
#cor(ds2[,1], ds2[,4])
#cor(ds2[,2], ds2[,3])

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

#Fig 5:

ggplot(data=dfPlot, aes(x=x, y=y, fill=Fit)) +
  geom_bar(stat="identity",position=position_dodge()) +
  coord_cartesian(ylim=c(0.8, 0.88)) +
  labs(title="Mean Effect of Regressing out Technical Covariates", y="Correlation, 10x vs bulk, log transformed data", x="Covariates regressed out") #+
  #theme(axis.text.x = element_text(angle=65, vjust=0.6))
  




  

#Fig 4: Create a combined plot of all covariates:
###########################################
plots = list(resCort1RemUMIFrac[[2]], resCort1GeneLength[[2]], resCort1GCFullLength[[2]], resCort1GCTail[[2]])
library("ggpubr")
dev.off()
gc()
fig = ggarrange( #when exporting this, make the y size larger(y=620)
  plotlist = plots, 
  nrow = 2, 
  ncol = 2,
  labels = c("A","B","C","D") 
) 
#check that the title is shown on the graph, it sometimes randomly disappears. 
annotate_figure(fig,
                top = text_grob("Log2 Fold Change Between 10x and Bulk vs Technical Covariates", face = "bold", size = 14))



