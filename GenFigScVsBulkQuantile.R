#Generates Fig. S3 in the paper. 
#Run GenerateScVsBulkDataQuantile.R before running this file!
#Note that this is a copy of the file GenFigScVsBulk.R, but here we use quantile normalization instead of TMM to 
#show that the UMICF regression cannot be explained by differences in gene expression histogram between single-cell and bulk.

library("ggplot2")
library("Seurat")

dataFolder = "C:/Work/R/RNASeqCTProfileEval/"
source(paste0(dataFolder, "FigureHelpFunc.R"))

#Gets Pearson correlation betweem UMI and Bulk for the EVAL dataset.
corrUMIVsBulkQ <- function(ds) {
  cor(ds$logUMIQ, ds$logBulkQ)
}

#Will regress out a fit in log space and update the fields "logUMIQ" and "LogUMIDivBulk"
#may produce some warnings about missing values that can be ignored
regrOutUMIVsBulkQ <- function(ds, fit) {
  ds2 = ds
  pred = predict(fit,ds2)
  fitNAFilter = !is.na(pred);
  #Now, regress out
  ds2$logUMIQ[fitNAFilter] = ds2$logUMIQ[fitNAFilter] - pred[fitNAFilter] + mean(pred[fitNAFilter])

  #also update the "LogUMIDivBulk"
  ds2$LogUMIDivBulk = ds2$logUMIQ - ds2$logBulkQ;
  
  return (ds2)
}

#Help function used from GenFigScVsBulk.
#Generates the plots in Fig 5 and some additional plots. The function operates on one covariate
#described by the params. Also returns values describing the improvement when regressing out the
#covariate.
genScToBulkCovGraphsQ <- function(ds, formulaUMI, formulaLFC, covIndex, covName, filter = NA) {
  if (!is.na(filter[1])) {
    ds = ds[filter,]
    print("filtering")
  }
  ind = sort(ds[,covIndex], index.return=T, na.last = T)
  dsSort = ds[ind$ix,];
  naFilt = !is.na(dsSort[,covIndex])
  numGenes = dim(dsSort)[1]
  plotFilter = !is.na(dsSort[,covIndex])
  dsPlot = data.frame(dsSort[plotFilter,covIndex], dsSort$logUMIQ[plotFilter])
  colnames(dsPlot) = c("x","y")
  loess_fit <- loess(formulaUMI, dsSort, span = 0.3)
  #sometimes have some NA, so we need to handle that
  dsLoess = data.frame(dsSort[naFilt,covIndex], predict(loess_fit, dsSort[naFilt,])) 
  colnames(dsLoess) = c("x","y")
  
  #Plot log expression in UMI data vs covariate 
  p1 = ggplot(dsPlot,aes(x=x,y=y))
  p1 = p1 + geom_point(alpha=0.3, shape=1)
  p1 = p1 + geom_line(data = dsLoess, colour="#FF0000", size=1.4, alpha=1)
  p1 = p1 + labs(y="10X gene expression (log2(pseudo-CPM))", x=covName)
  #p1<-p1 + labs(title="Visualization of Batch Effects")
  #p<-p + coord_cartesian(xlim=c(-0.09, 1), ylim=c(-0.09, 1))#create room for PC label
  #p<-p + theme( axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) 
  print(p1)
  
  #Plot log fold change between UMI and bulk vs covariate 
  #plot(dsSort[,covIndex], dsSort$LogUMIDivBulk)
  loess_fit2 = loess(formulaLFC, dsSort)
  #lines(dsSort[,covIndex], predict(loess_fit2), col = "red")
  #make linear fit as well
  lm2 = lm(formulaLFC, dsSort)
  #lines(dsSort[,covIndex], predict(lm2), col = "blue")
  
  #Plot log fold change between UMI and bulk vs covariate 
  dsPlot = data.frame(dsSort[plotFilter,covIndex], dsSort$LogUMIDivBulk[plotFilter])
  colnames(dsPlot) = c("x","y")
  
  dsLoess = data.frame(dsSort[naFilt,covIndex], predict(loess_fit2, dsSort[naFilt,]))
  colnames(dsLoess) = c("x","y")
  dsLin = data.frame(dsSort[naFilt,covIndex], predict(lm2, dsSort[naFilt,]))
  colnames(dsLin) = c("x","y")
  
  p2 = ggplot(dsPlot,aes(x=x,y=y))
  p2 = p2 + geom_point(alpha=0.3, shape=1)
  p2 = p2 + geom_line(data = dsLoess, colour="#FF0000", alpha=1, size=1.4)
  p2 = p2 + geom_line(data = dsLin, colour="#00BB00", alpha=1, size=1.4)
  p2 = p2 + labs(y="Log2 fold change, 10X vs bulk", x=covName)
  #p1<-p1 + labs(title="Visualization of Batch Effects")
  #p<-p + coord_cartesian(xlim=c(-0.09, 1), ylim=c(-0.09, 1))#create room for PC label
  #p<-p + theme( axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) 
  print(p2)
  
  
  #Regress out covariate:
  dsRegr = regrOutUMIVsBulkQ(dsSort, loess_fit2)
  dsRegrLin = regrOutUMIVsBulkQ(dsSort, lm2)
  uvsb = corrUMIVsBulkQ(dsSort)
  print(uvsb)
  ruvsb = corrUMIVsBulkQ(dsRegr)
  print(ruvsb)
  ruvsbl = corrUMIVsBulkQ(dsRegrLin)
  print(ruvsbl)
  #plot(dsRegr[,covIndex], dsRegr$LogUMIDivBulk) # for test only
  
  res = list(p1, p2, uvsb, ruvsb, ruvsbl)
  
  return (res)
}




###########################
## Technical biases between single-cell and bulk
###########################

#Four things: 1. UMIs vs counts (removed counts' fraction), 2. gene length, 3. GC content 4. GC content tail

cort1 = read.table(paste0(dataFolder, "/data/ScVsBulkQuantileCortex1.txt"), sep="\t")
cort2 = read.table(paste0(dataFolder, "/data/ScVsBulkQuantileCortex2.txt"), sep="\t")

#check correlation between UMIFrac for cort1 and cort2. Not super, but it seems it is there
plot(cort1$remUMIFrac, cort2$remUMIFrac)

#Remove all lowly expressed genes, so much noise there
#filt = (cort1$logUMIQ >= log2(1.05)) & (cort1$logBulkQ >= log2(1.05))   #filter on pseudo-TPM of 1
#filtCortExpr1 = cort1[filt,]
#filt2 = (cort2$logUMIQ >= log2(1.05)) & (cort2$logBulkQ >= log2(1.05)) #filter on pseudo-TPM of 1
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

sum(is.na(filtCort1$remUMIFrac))
sum(is.na(filtCort2$remUMIFrac))


#experiment with not using umi copy fraction for genes with very few counts - it seems 5 UMIs or less gives only noise
filtCortRemRUFExpr1 = filtCortExpr1;
filtCortRemRUFExpr1$remUMIFrac[filtCortExpr2$UMI < 6] = NA
#filtCortRemRUFExpr1$remUMIFrac[filtCortExpr1$UMITPM < 0.001 | filtCortExpr1$bulkTPM < 0.001] = NA
filtCortRemRUF1 = filtCortRemRUFExpr1[filtAll1,]
#resCort1RemUMIFracSpec = genScToBulkCovGraphsQ(filtCortRemRUF1, logUMIQ ~ remUMIFrac, LogUMIDivBulk ~ remUMIFrac, 7, "UMI copy fraction")

filtCortRemRUFExpr2 = filtCortExpr2;
filtCortRemRUFExpr2$remUMIFrac[filtCortExpr1$UMI < 6] = NA
#filtCortRemRUFExpr1$remUMIFrac[filtCortExpr1$UMITPM < 0.001 | filtCortExpr1$bulkTPM < 0.001] = NA
filtCortRemRUF2 = filtCortRemRUFExpr2[filtAll2,]
#resCort2RemUMIFracSpec = genScToBulkCovGraphsQ(filtCortRemRUF2, logUMIQ ~ remUMIFrac, LogUMIDivBulk ~ remUMIFrac, 7, "UMI copy fraction")


#0 Copies per UMI
###########################
#resCort1CopiesPerUMI = genScToBulkCovGraphsQ(filtCort1, logUMIQ ~ CopiesPerUMIOtherSample, LogUMIDivBulk ~ CopiesPerUMIOtherSample, 14, "Copies per UMI")
#resCort2RemUMIFrac = genScToBulkCovGraphsQ(filtCort2, logUMIQ ~ remUMIFrac, LogUMIDivBulk ~ remUMIFrac, 7, "UMI copy fraction")


#1 Removed counts' fraction
###########################
#resCort1RemUMIFrac = genScToBulkCovGraphsQ(filtCort1, logUMIQ ~ remUMIFrac, LogUMIDivBulk ~ remUMIFrac, 7, "UMI copy fraction")
#resCort2RemUMIFrac = genScToBulkCovGraphsQ(filtCort2, logUMIQ ~ remUMIFrac, LogUMIDivBulk ~ remUMIFrac, 7, "UMI copy fraction")
#use the data where UMI copy fraction is removed where it is based on too few molecules
resCort1RemUMIFrac = genScToBulkCovGraphsQ(filtCortRemRUF1, logUMIQ ~ remUMIFrac, LogUMIDivBulk ~ remUMIFrac, 6, "UMI copy fraction")
resCort2RemUMIFrac = genScToBulkCovGraphsQ(filtCortRemRUF2, logUMIQ ~ remUMIFrac, LogUMIDivBulk ~ remUMIFrac, 6, "UMI copy fraction")

resCort1RemUMIFrac2 = genScToBulkCovGraphsQ(filtCortRemRUF1, logUMIQ ~ remUMIFrac, LogUMIDivBulk ~ remUMIFrac, 6, "UMI copy fraction")


#2. Gene length
###########################

#filter out outlier genes above 10000 in length
resCort1GeneLength = genScToBulkCovGraphsQ(filtCortRemRUF1, logUMIQ ~ geneLength, LogUMIDivBulk ~ geneLength, 5, "Transcript length")
resCort2GeneLength = genScToBulkCovGraphsQ(filtCortRemRUF2, logUMIQ ~ geneLength, LogUMIDivBulk ~ geneLength, 5, "Transcript length")

#3. GC Content full length
###########################
resCort1GCFullLength = genScToBulkCovGraphsQ(filtCortRemRUF1, logUMIQ ~ gcFullLength, LogUMIDivBulk ~ gcFullLength, 3, "GC content entire transcript")
resCort2GCFullLength = genScToBulkCovGraphsQ(filtCortRemRUF2, logUMIQ ~ gcFullLength, LogUMIDivBulk ~ gcFullLength, 3, "GC content entire transcript")

#4. GC Content tail
###########################
resCort1GCTail = genScToBulkCovGraphsQ(filtCortRemRUF1, logUMIQ ~ gcTail, LogUMIDivBulk ~ gcTail, 4, "GC content transcript tail")
resCort2GCTail = genScToBulkCovGraphsQ(filtCortRemRUF2, logUMIQ ~ gcTail, LogUMIDivBulk ~ gcTail, 4, "GC content transcript tail")

#5. Test to regress out all covariates:
###########################


#Regress out all covariates
formAll = LogUMIDivBulk ~ remUMIFrac + geneLength + gcFullLength + gcTail
loess_fitAll1 <- loess(formAll, filtCortRemRUF1)
lmAll1 <- lm(formAll, filtCortRemRUF1)
loess_fitAll2 <- loess(formAll, filtCortRemRUF2)
lmAll2 <- lm(formAll, filtCortRemRUF2)

cort1RegrAllLoess = regrOutUMIVsBulkQ(filtCortRemRUF1, loess_fitAll1)
cort1RegrAllLin = regrOutUMIVsBulkQ(filtCortRemRUF1, lmAll1)
cort2RegrAllLoess = regrOutUMIVsBulkQ(filtCortRemRUF2, loess_fitAll2)
cort2RegrAllLin = regrOutUMIVsBulkQ(filtCortRemRUF2, lmAll2)

corAllLoess1 = corrUMIVsBulkQ(cort1RegrAllLoess)
corAllLin1 = corrUMIVsBulkQ(cort1RegrAllLin)
corAllLoess2 = corrUMIVsBulkQ(cort2RegrAllLoess)
corAllLin2 = corrUMIVsBulkQ(cort2RegrAllLin)

#regress out all but gcTail
formAllButGCTail = LogUMIDivBulk ~ remUMIFrac + geneLength + gcFullLength
loess_fitAllButGCTail1 <- loess(formAllButGCTail, filtCortRemRUF1)
lmAllButGCTail1 <- lm(formAllButGCTail, filtCortRemRUF1)
loess_fitAllButGCTail2 <- loess(formAllButGCTail, filtCortRemRUF2)
lmAllButGCTail2 <- lm(formAllButGCTail, filtCortRemRUF2)

cort1RegrAllButGCTailLoess = regrOutUMIVsBulkQ(filtCortRemRUF1, loess_fitAllButGCTail1)
cort1RegrAllButGCTailLin = regrOutUMIVsBulkQ(filtCortRemRUF1, lmAllButGCTail1)
cort2RegrAllButGCTailLoess = regrOutUMIVsBulkQ(filtCortRemRUF2, loess_fitAllButGCTail2)
cort2RegrAllButGCTailLin = regrOutUMIVsBulkQ(filtCortRemRUF2, lmAllButGCTail2)

corAllButGCTailLoess1 = corrUMIVsBulkQ(cort1RegrAllButGCTailLoess)
corAllButGCTailLin1 = corrUMIVsBulkQ(cort1RegrAllButGCTailLin)
corAllButGCTailLoess2 = corrUMIVsBulkQ(cort2RegrAllButGCTailLoess)
corAllButGCTailLin2 = corrUMIVsBulkQ(cort2RegrAllButGCTailLin)

#regress out remUMIFrac and gcFullLength (i.e. remove gene length as well)
formRemUMIFracAndGCFullLength = LogUMIDivBulk ~ remUMIFrac + gcFullLength
loess_fitRemUMIFracAndGCFullLength1 <- loess(formRemUMIFracAndGCFullLength, filtCortRemRUF1)
lmRemUMIFracAndGCFullLength1 <- lm(formRemUMIFracAndGCFullLength, filtCortRemRUF1)
loess_fitRemUMIFracAndGCFullLength2 <- loess(formRemUMIFracAndGCFullLength, filtCortRemRUF2)
lmRemUMIFracAndGCFullLength2 <- lm(formRemUMIFracAndGCFullLength, filtCortRemRUF2)

cort1RegrRemUMIFracAndGCFullLengthLoess = regrOutUMIVsBulkQ(filtCortRemRUF1, loess_fitRemUMIFracAndGCFullLength1)
cort1RegrRemUMIFracAndGCFullLengthLin = regrOutUMIVsBulkQ(filtCortRemRUF1, lmRemUMIFracAndGCFullLength1)
cort2RegrRemUMIFracAndGCFullLengthLoess = regrOutUMIVsBulkQ(filtCortRemRUF2, loess_fitRemUMIFracAndGCFullLength2)
cort2RegrRemUMIFracAndGCFullLengthLin = regrOutUMIVsBulkQ(filtCortRemRUF2, lmRemUMIFracAndGCFullLength2)

corRemUMIFracAndGCFullLengthLoess1 = corrUMIVsBulkQ(cort1RegrRemUMIFracAndGCFullLengthLoess)
corRemUMIFracAndGCFullLengthLin1 = corrUMIVsBulkQ(cort1RegrRemUMIFracAndGCFullLengthLin)
corRemUMIFracAndGCFullLengthLoess2 = corrUMIVsBulkQ(cort2RegrRemUMIFracAndGCFullLengthLoess)
corRemUMIFracAndGCFullLengthLin2 = corrUMIVsBulkQ(cort2RegrRemUMIFracAndGCFullLengthLin)

formGeneLengthAndGCFullLength = LogUMIDivBulk ~ geneLength + gcFullLength
loess_fitGeneLengthAndGCFullLength1 <- loess(formGeneLengthAndGCFullLength, filtCortRemRUF1)
cort1GeneLengthAndGCFullLengthLoess = regrOutUMIVsBulkQ(filtCortRemRUF1, loess_fitGeneLengthAndGCFullLength1)
corGeneLengthAndGCFullLengthLoess1 = corrUMIVsBulkQ(cort1GeneLengthAndGCFullLengthLoess)

plot(cort1GeneLengthAndGCFullLengthLoess$logBulkQ, cort1GeneLengthAndGCFullLengthLoess$logUMIQ)
plot(filtCortRemRUF1$logBulkQ, filtCortRemRUF1$logUMIQ)

#check correlation between cortex 1 and 2 as comparison
print("Correlation values in the text");
ds2 = cbind(filtCortRemRUF1$logBulkQ, filtCortRemRUF1$logUMIQ, filtCortRemRUF2$logBulkQ, filtCortRemRUF2$logUMIQ)
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

x = factor(1:8, 1:8, c("None", "UMICF", "Tr. length", "GC cont.", "GC cont. tail", "UMICF\nGC cont.\nTr. length\nGC cont. tail", "UMICF\nGC cont.\nTr. length", "UMICF\nGC cont."))
Fit = factor(rep(c(1,2), each=length(meanLin)), 1:2, c("linear", "loess"))
y = c(meanLin, meanLoess)

dfPlot = data.frame(x, y, Fit)

#Fig S3 (Fig 6 for quantile):
###########################################

bp = ggplot(data=dfPlot, aes(x=x, y=y, fill=Fit)) +
  geom_bar(stat="identity",position=position_dodge()) +
  coord_cartesian(ylim=c(0.8, 0.88)) +
  labs( y="Correlation - 10x vs bulk", x="Covariates regressed out") +
  scale_fill_manual(values=c("#ACBEE8", "#6382D3")) +
  theme(axis.text = element_text(size = 8))
#+
  #theme(axis.text.x = element_text(angle=65, vjust=0.6))
  
p1 = plotCorr(filtCortRemRUF1$logUMIQ, filtCortRemRUF1$logBulkQ)
p2 = plotCorr(cort1RegrRemUMIFracAndGCFullLengthLoess$logUMIQ, cort1RegrRemUMIFracAndGCFullLengthLoess$logBulkQ)

library("ggpubr")


fig6QComp = ggarrange( #when exporting this, make the x size larger(x=800)
  ggarrange(p1,p2,ncol=2, labels=c("A","B")),
  bp, 
  nrow = 2, 
  ncol = 1,
  labels = c("", "C") 
) 

#check that the title is shown on the graph, it sometimes randomly disappears. 
figS3 = fig6QComp
#figS3 = annotate_figure(fig6QComp,
#                top = text_grob("Effect of Regressing out Technical Covariates - Quantile Norm.", face = "bold", size = 14))

ggsave(
  paste0(fig_path, "figS3.png"),
  plot = figS3, device = "png",
  width = 6, height = 6, dpi = 300)

#Not used (Fig 5 for quantile): Create a combined plot of all covariates:
###########################################
plots = list(resCort1RemUMIFrac[[2]], resCort1GeneLength[[2]], resCort1GCFullLength[[2]], resCort1GCTail[[2]])
dev.off()
gc()
fig = ggarrange( #when exporting this, make the x size larger(x=800)
  plotlist = plots, 
  nrow = 2, 
  ncol = 2,
  labels = c("A","B","C","D") 
) 
#check that the title is shown on the graph, it sometimes randomly disappears. 
annotate_figure(fig,
                top = text_grob("Log2 Fold Change Between 10x and Bulk vs Technical Covariates", face = "bold", size = 14))



#Not used : Fig S2: Gene expression vs UMI Copy Fraction
###########################################
figS2 = resCort1RemUMIFrac[[1]] + labs( title="Gene Expression vs UMI Copy Fraction")
figS2

#Some tests
############################################


# TC001B: Check that regressing out a covariate leads to that the correlation with the covariate is lost
form1 = LogUMIDivBulk ~ remUMIFrac
loess1 <- loess(form1, filtCortRemRUF1)
lm1 <- lm(form1, filtCortRemRUF1)
regrOut1Loess = regrOutUMIVsBulkQ(filtCortRemRUF1, loess1)
regrOut1Lin = regrOutUMIVsBulkQ(filtCortRemRUF1, lm1)

plot(filtCortRemRUF1$remUMIFrac, filtCortRemRUF1$LogUMIDivBulk)
abline(lm1,col="red")

plot(regrOut1Loess$remUMIFrac, regrOut1Loess$LogUMIDivBulk)
lm2 <- lm(form1, regrOut1Loess)
abline(lm2,col="green") # should be reasonably flat

plot(regrOut1Lin$remUMIFrac, regrOut1Lin$LogUMIDivBulk)
lm3 <- lm(form1, regrOut1Lin)
abline(lm3,col="blue") # should be completely flat

#All looks fine!

