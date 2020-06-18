#Generates Fig. 5 and 6 in the paper. 
#Run GenerateScVsBulkData.R before running this file!

library("ggplot2")
library("Seurat")

dataFolder = "C:/Work/R/RNASeqCTProfileEval/"
source(paste0(dataFolder, "FigureHelpFunc.R"))

fig_path = "Z:/projects/Cell type profiles/figures/"


#Gets Pearson correlation betweem UMI and Bulk for the EVAL dataset.
corrUMIVsBulk <- function(ds) {
  cor(ds$logUMITMM, ds$logBulkTMM)
}

#so, this function creates 10000 bootstraps and calculates the correlation each time
corrUMIVsBulkBootstrap <- function(ds, bootstraps) {
  numGenes = dim(ds)[1]
  corVals = rep(0, 10000);
  for (i in 1:10000) {
    sel = bootstraps[,i] #sample(numGenes, size=numGenes, replace=T)
    corVals[i] = cor(ds$logUMITMM[sel], ds$logBulkTMM[sel])
  }
  return (corVals)
}

getUpperConfInterval <- function(bootstrapData) {
  return (sort(bootstrapData, decreasing=T)[round(length(bootstrapData)*0.025)])
}
getLowerConfInterval <- function(bootstrapData) {
  return (sort(bootstrapData, decreasing=F)[round(length(bootstrapData)*0.025)])
}

#test the bootstrapping
#a = rnorm(numGenes, mean = 100, sd = 10)
#b = rnorm(numGenes, mean = 100, sd = 10)
#plot(a,b)
#cor(a,b)
#dsTest = data.frame(logUMITMM=a, logBulkTMM=b)
#cors = corrUMIVsBulkBootstrap(dsTest, bootstraps)
#hist(cors)
#getUpperConfInterval(cors)
#getLowerConfInterval(cors)#these look reasonable


#Will regress out a fit in log space and update the fields "logUMITMM" and "LogUMIDivBulk"
#may produce some warnings about missing values that can be ignored
regrOutUMIVsBulk <- function(ds, fit) {
  ds2 = ds
  pred = predict(fit,ds2)
  fitNAFilter = !is.na(pred);
  #Now, regress out
  ds2$logUMITMM[fitNAFilter] = ds2$logUMITMM[fitNAFilter] - pred[fitNAFilter] + mean(pred[fitNAFilter])
  
  #also update the "LogUMIDivBulk"
  ds2$LogUMIDivBulk = ds2$logUMITMM - ds2$logBulkTMM;
  
  return (ds2)
}


#Help function used from GenFigScVsBulk.
#Generates the plots in Fig 5 and some additional plots. The function operates on one covariate
#described by the params. Also returns values describing the improvement when regressing out the
#covariate.
genScToBulkCovGraphs <- function(ds, formulaUMI, formulaLFC, covIndex, covName, bootstraps, filter = NA) {
  if (!is.na(filter[1])) {
    ds = ds[filter,]
    print("filtering")
  }
  ind = sort(ds[,covIndex], index.return=T, na.last = T)
  dsSort = ds[ind$ix,];
  naFilt = !is.na(dsSort[,covIndex])
  numGenes = dim(dsSort)[1]
  plotFilter = !is.na(dsSort[,covIndex])
  dsPlot = data.frame(dsSort[plotFilter,covIndex], dsSort$logUMITMM[plotFilter])
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
  #plot(dsSort[,covIndex], dsSort$logUMITMM)
  #loess_fit <- loess(formulaUMI, dsSort)
  #lines(dsSort[,covIndex], predict(loess_fit), col = "red")
  
  #Plot log fold change between UMI and bulk vs covariate 
  #plot(dsSort[,covIndex], dsSort$LogUMIDivBulk)
  loess_fit2 = loess(formulaLFC, dsSort)
  #lines(dsSort[,covIndex], predict(loess_fit2), col = "red")
  #make linear fit as well
  lm2 = lm(formulaLFC, dsSort)
  print(lm2)
  print(summary(lm2))
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
  dsRegr = regrOutUMIVsBulk(dsSort, loess_fit2)
  dsRegrLin = regrOutUMIVsBulk(dsSort, lm2)
  uvsb = corrUMIVsBulk(dsSort)
  print(uvsb)
  ruvsb = corrUMIVsBulk(dsRegr)
  ruvsbbootstrap = corrUMIVsBulkBootstrap(dsRegr, bootstraps)
  print(ruvsb)
  ruvsbl = corrUMIVsBulk(dsRegrLin)
  ruvsblbootstrap = corrUMIVsBulkBootstrap(dsRegrLin, bootstraps)
  print(ruvsbl)
  #plot(dsRegr[,covIndex], dsRegr$LogUMIDivBulk) # for test only
  
  res = list(p1, p2, uvsb, ruvsb, ruvsbl, ruvsbbootstrap, ruvsblbootstrap)#the two last are pairs of confidence intervals
  
  return (res)
}




###########################
## Figure 5 and 6 - Technical biases between single-cell and bulk
###########################

#Four things: 1. UMIs vs counts (removed counts' fraction), 2. gene length, 3. GC content 4. GC content tail

cort1 = read.table(paste0(dataFolder, "/data/ScVsBulkCortex1.txt"), sep="\t")
cort2 = read.table(paste0(dataFolder, "/data/ScVsBulkCortex2.txt"), sep="\t")

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

sum(is.na(filtCort1$remUMIFrac))
sum(is.na(filtCort2$remUMIFrac))


#experiment with not using umi copy fraction for genes with very few counts - it seems 5 UMIs or less gives only noise
filtCortRemRUFExpr1 = filtCortExpr1;
filtCortRemRUFExpr1$remUMIFrac[filtCortExpr1$UMI < 6] = NA
filtCortRemRUF1 = filtCortRemRUFExpr1[filtAll1,]

filtCortRemRUFExpr2 = filtCortExpr2;
filtCortRemRUFExpr2$remUMIFrac[filtCortExpr2$UMI < 6] = NA
filtCortRemRUF2 = filtCortRemRUFExpr2[filtAll2,]

numGenes = length(filtCortRemRUF1$UMI) #same for both samples

#Generate bootstraps for correlation
bootstraps = matrix(0, nrow = numGenes, ncol = 10000);
for (i in 1:10000) {
  bootstraps[,i] = sample(numGenes, size=numGenes, replace=T)
}




#1 Removed counts' fraction
###########################
#use the data where UMI copy fraction is removed where it is based on too few molecules
resCort1RemUMIFrac = genScToBulkCovGraphs(filtCortRemRUF1, logUMITMM ~ remUMIFrac, LogUMIDivBulk ~ remUMIFrac, 6, "UMI copy fraction", bootstraps)
resCort2RemUMIFrac = genScToBulkCovGraphs(filtCortRemRUF2, logUMITMM ~ remUMIFrac, LogUMIDivBulk ~ remUMIFrac, 6, "UMI copy fraction", bootstraps)


#2. Gene length
###########################

#filter out outlier genes above 10000 in length
resCort1GeneLength = genScToBulkCovGraphs(filtCortRemRUF1, logUMITMM ~ geneLength, LogUMIDivBulk ~ geneLength, 5, "Transcript length", bootstraps)
resCort2GeneLength = genScToBulkCovGraphs(filtCortRemRUF2, logUMITMM ~ geneLength, LogUMIDivBulk ~ geneLength, 5, "Transcript length", bootstraps)

#3. GC Content full length
###########################
resCort1GCFullLength = genScToBulkCovGraphs(filtCortRemRUF1, logUMITMM ~ gcFullLength, LogUMIDivBulk ~ gcFullLength, 3, "GC content entire transcript", bootstraps)
resCort2GCFullLength = genScToBulkCovGraphs(filtCortRemRUF2, logUMITMM ~ gcFullLength, LogUMIDivBulk ~ gcFullLength, 3, "GC content entire transcript", bootstraps)

#4. GC Content tail
###########################
resCort1GCTail = genScToBulkCovGraphs(filtCortRemRUF1, logUMITMM ~ gcTail, LogUMIDivBulk ~ gcTail, 4, "GC content transcript tail", bootstraps)
resCort2GCTail = genScToBulkCovGraphs(filtCortRemRUF2, logUMITMM ~ gcTail, LogUMIDivBulk ~ gcTail, 4, "GC content transcript tail", bootstraps)

#5. Test to regress out all covariates:
###########################


#Regress out all covariates
formAll = LogUMIDivBulk ~ remUMIFrac + geneLength + gcFullLength + gcTail
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

bootstrAllLoess1 = corrUMIVsBulkBootstrap(cort1RegrAllLoess, bootstraps)
bootstrAllLin1 = corrUMIVsBulkBootstrap(cort1RegrAllLin, bootstraps)
bootstrAllLoess2 = corrUMIVsBulkBootstrap(cort2RegrAllLoess, bootstraps)
bootstrAllLin2 = corrUMIVsBulkBootstrap(cort2RegrAllLin, bootstraps)



#regress out all but gcTail
formAllButGCTail = LogUMIDivBulk ~ remUMIFrac + geneLength + gcFullLength
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

bootstrAllButGCTailLoess1 = corrUMIVsBulkBootstrap(cort1RegrAllButGCTailLoess, bootstraps)
bootstrAllButGCTailLin1 = corrUMIVsBulkBootstrap(cort1RegrAllButGCTailLin, bootstraps)
bootstrAllButGCTailLoess2 = corrUMIVsBulkBootstrap(cort2RegrAllButGCTailLoess, bootstraps)
bootstrAllButGCTailLin2 = corrUMIVsBulkBootstrap(cort2RegrAllButGCTailLin, bootstraps)


#regress out remUMIFrac and gcFullLength (i.e. remove gene length as well)
formRemUMIFracAndGCFullLength = LogUMIDivBulk ~ remUMIFrac + gcFullLength
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

bootstrRemUMIFracAndGCFullLengthLoess1 = corrUMIVsBulkBootstrap(cort1RegrRemUMIFracAndGCFullLengthLoess, bootstraps)
bootstrRemUMIFracAndGCFullLengthLin1 = corrUMIVsBulkBootstrap(cort1RegrRemUMIFracAndGCFullLengthLin, bootstraps)
bootstrRemUMIFracAndGCFullLengthLoess2 = corrUMIVsBulkBootstrap(cort2RegrRemUMIFracAndGCFullLengthLoess, bootstraps)
bootstrRemUMIFracAndGCFullLengthLin2 = corrUMIVsBulkBootstrap(cort2RegrRemUMIFracAndGCFullLengthLin, bootstraps)


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
#meanLin = (c1Lin + c2Lin)/2
c1Loess = c(resCort1RemUMIFrac[[3]], resCort1RemUMIFrac[[4]], resCort1GeneLength[[4]], resCort1GCFullLength[[4]],
            resCort1GCTail[[4]], corAllLoess1, corAllButGCTailLoess1, corRemUMIFracAndGCFullLengthLoess1)
c2Loess = c(resCort2RemUMIFrac[[3]], resCort2RemUMIFrac[[4]], resCort2GeneLength[[4]], resCort2GCFullLength[[4]],
            resCort2GCTail[[4]], corAllLoess2, corAllButGCTailLoess2, corRemUMIFracAndGCFullLengthLoess2)
#meanLoess = (c1Loess + c2Loess)/2
allData = c(c1Lin, c2Lin, c1Loess, c2Loess)


corTable = data.frame(c1Lin, c2Lin, meanLin, c1Loess, c2Loess, meanLoess)
round(corTable,3)

x = factor(rep(1:8,4), 1:8, c("None", "UMICF", "Tr. length", "GC cont.", "GC cont. tail", "UMICF\nGC cont.\nTr. length\nGC cont. tail", "UMICF\nGC cont.\nTr. length", "UMICF\nGC cont."))
Fit = factor(rep(c(1,2,3,4), each=length(c1Lin)), 1:4, c("lin. s1", "lin. s2", "loess s1", "loess s2"))
#y = c(meanLin, meanLoess)
y = allData

#create bootstrapped uncertainties
#get bootstrap correlations for the "none case"
bootstrNone1 = corrUMIVsBulkBootstrap(filtCortRemRUF1, bootstraps)
bootstrNone2 = corrUMIVsBulkBootstrap(filtCortRemRUF2, bootstraps)

ub1Lin = c(getUpperConfInterval(bootstrNone1), getUpperConfInterval(resCort1RemUMIFrac[[7]]),
           getUpperConfInterval(resCort1GeneLength[[7]]), getUpperConfInterval(resCort1GCFullLength[[7]]), 
           getUpperConfInterval(resCort1GCTail[[7]]), getUpperConfInterval(bootstrAllLin1), 
           getUpperConfInterval(bootstrAllButGCTailLin1), getUpperConfInterval(bootstrRemUMIFracAndGCFullLengthLin1) )
ub1Loess = c(getUpperConfInterval(bootstrNone1), getUpperConfInterval(resCort1RemUMIFrac[[6]]),
           getUpperConfInterval(resCort1GeneLength[[6]]), getUpperConfInterval(resCort1GCFullLength[[6]]), 
           getUpperConfInterval(resCort1GCTail[[6]]), getUpperConfInterval(bootstrAllLoess1), 
           getUpperConfInterval(bootstrAllButGCTailLoess1), getUpperConfInterval(bootstrRemUMIFracAndGCFullLengthLoess1) )
ub2Lin = c(getUpperConfInterval(bootstrNone2), getUpperConfInterval(resCort2RemUMIFrac[[7]]),
           getUpperConfInterval(resCort2GeneLength[[7]]), getUpperConfInterval(resCort2GCFullLength[[7]]), 
           getUpperConfInterval(resCort2GCTail[[7]]), getUpperConfInterval(bootstrAllLin2), 
           getUpperConfInterval(bootstrAllButGCTailLin2), getUpperConfInterval(bootstrRemUMIFracAndGCFullLengthLin2) )
ub2Loess = c(getUpperConfInterval(bootstrNone2), getUpperConfInterval(resCort2RemUMIFrac[[6]]),
             getUpperConfInterval(resCort2GeneLength[[6]]), getUpperConfInterval(resCort2GCFullLength[[6]]), 
             getUpperConfInterval(resCort2GCTail[[6]]), getUpperConfInterval(bootstrAllLoess2), 
             getUpperConfInterval(bootstrAllButGCTailLoess2), getUpperConfInterval(bootstrRemUMIFracAndGCFullLengthLoess2) )
lb1Lin = c(getLowerConfInterval(bootstrNone1), getLowerConfInterval(resCort1RemUMIFrac[[7]]),
           getLowerConfInterval(resCort1GeneLength[[7]]), getLowerConfInterval(resCort1GCFullLength[[7]]), 
           getLowerConfInterval(resCort1GCTail[[7]]), getLowerConfInterval(bootstrAllLin1), 
           getLowerConfInterval(bootstrAllButGCTailLin1), getLowerConfInterval(bootstrRemUMIFracAndGCFullLengthLin1) )
lb1Loess = c(getLowerConfInterval(bootstrNone1), getLowerConfInterval(resCort1RemUMIFrac[[6]]),
             getLowerConfInterval(resCort1GeneLength[[6]]), getLowerConfInterval(resCort1GCFullLength[[6]]), 
             getLowerConfInterval(resCort1GCTail[[6]]), getLowerConfInterval(bootstrAllLoess1), 
             getLowerConfInterval(bootstrAllButGCTailLoess1), getLowerConfInterval(bootstrRemUMIFracAndGCFullLengthLoess1) )
lb2Lin = c(getLowerConfInterval(bootstrNone2), getLowerConfInterval(resCort2RemUMIFrac[[7]]),
           getLowerConfInterval(resCort2GeneLength[[7]]), getLowerConfInterval(resCort2GCFullLength[[7]]), 
           getLowerConfInterval(resCort2GCTail[[7]]), getLowerConfInterval(bootstrAllLin2), 
           getLowerConfInterval(bootstrAllButGCTailLin2), getLowerConfInterval(bootstrRemUMIFracAndGCFullLengthLin2) )
lb2Loess = c(getLowerConfInterval(bootstrNone2), getLowerConfInterval(resCort2RemUMIFrac[[6]]),
             getLowerConfInterval(resCort2GeneLength[[6]]), getLowerConfInterval(resCort2GCFullLength[[6]]), 
             getLowerConfInterval(resCort2GCTail[[6]]), getLowerConfInterval(bootstrAllLoess2), 
             getLowerConfInterval(bootstrAllButGCTailLoess2), getLowerConfInterval(bootstrRemUMIFracAndGCFullLengthLoess2) )
ub = c(ub1Lin, ub2Lin, ub1Loess, ub2Loess)
lb = c(lb1Lin, lb2Lin, lb1Loess, lb2Loess)


dfPlot = data.frame(x, y, Fit, lb, ub)

#Fig 6:
###########################################

bp = ggplot(data=dfPlot, aes(x=x, y=y, fill=Fit)) +
  geom_bar(stat="identity",position=position_dodge()) +
  geom_errorbar(aes(ymin=lb, ymax=ub), width=.2,
                position=position_dodge(.9)) +
  coord_cartesian(ylim=c(0.8, 0.89)) +
  labs( y="Correlation - 10x vs bulk", x="Covariates regressed out") +
  #scale_fill_manual(values=c("#ACBEE8", "#6382D3", "#82E182", "#229A22")) +
  scale_fill_manual(values=c("#82E182", "#229A22", "#E47060", "#AC210E")) + #so, use red and green as in figure 5
  theme(axis.text = element_text(size = 8))
#+
  #theme(axis.text.x = element_text(angle=65, vjust=0.6))

p1 = plotCorr(filtCortRemRUF1$logUMITMM, filtCortRemRUF1$logBulkTMM)
p2 = plotCorr(cort1RegrRemUMIFracAndGCFullLengthLoess$logUMITMM, cort1RegrRemUMIFracAndGCFullLengthLoess$logBulkTMM)

library("ggpubr")


fig6Comb = ggarrange( #when exporting this, make the x size larger(x=800)
  ggarrange(p1,p2,ncol=2, labels=c("A","B")),
  bp, 
  nrow = 2, 
  ncol = 1,
  labels = c("", "C") 
) 

#check that the title is shown on the graph, it sometimes randomly disappears. 
fig6 = fig6Comb
#fig6 = annotate_figure(fig6Comb,
#                top = text_grob("Effect of Regressing out Technical Covariates", face = "bold", size = 14))

fig6

ggsave(
  paste0(fig_path, "Fig6.png"),
  plot = fig6, device = "png",
  width = 6, height = 6, dpi = 300)


#Fig 5: Create a combined plot of all covariates:
###########################################
plots = list(resCort1RemUMIFrac[[2]], resCort1GeneLength[[2]], resCort1GCFullLength[[2]], resCort1GCTail[[2]])
dev.off()
gc()
fig5Comb = ggarrange( #when exporting this, make the x size larger(x=800)
  plotlist = plots, 
  nrow = 2, 
  ncol = 2,
  labels = c("A","B","C","D") 
) 
#check that the title is shown on the graph, it sometimes randomly disappears. 
fig5 = fig5Comb #skip title
#fig5 = annotate_figure(fig5Comb,
#                top = text_grob("Log2 Fold Change Between 10x and Bulk vs Technical Covariates", face = "bold", size = 14))

fig5

ggsave(
  paste0(fig_path, "Fig5.png"),
  plot = fig5, device = "png",
  width = 6, height = 6, dpi = 300)

#so, the bootstrapped covariates are generated with the same gene bootstrap matrix, so the comparison
#shall be paired (just bootstrapping may change the correlation depending on gene selection, but this
#way that will be the same for both pairs being compared)

#Wilcoxon signed rank tests that UMICF is a better covariate than GC content
wilc1lin = wilcox.test(x = resCort1RemUMIFrac[[7]], y = resCort1GCFullLength[[7]],
                    alternative = "greater", paired = T)
wilc1loess = wilcox.test(x = resCort1RemUMIFrac[[6]], y = resCort1GCFullLength[[6]],
                    alternative = "greater", paired = T)
wilc2lin = wilcox.test(x = resCort2RemUMIFrac[[7]], y = resCort2GCFullLength[[7]],
                    alternative = "greater", paired = T)
wilc2loess = wilcox.test(x = resCort2RemUMIFrac[[6]], y = resCort2GCFullLength[[6]],
                    alternative = "greater", paired = T)


#Fig S2: Gene expression vs UMI Copy Fraction
###########################################
figS2a = resCort1RemUMIFrac[[1]]
figS2 = figS2a # skip title
#figS2 = annotate_figure(figS2a,
#                        top = text_grob("Gene Expression vs UMI Copy Fraction", face = "bold", size = 14))
figS2
ggsave(
  paste0(fig_path, "FigS2.png"),
  plot = figS2, device = "png",
  width = 6, height = 6, dpi = 300)

#Some tests
############################################


# TC001: Check that regressing out a covariate leads to that the correlation with the covariate is lost
form1 = LogUMIDivBulk ~ remUMIFrac
loess1 <- loess(form1, filtCortRemRUF1)
lm1 <- lm(form1, filtCortRemRUF1)
regrOut1Loess = regrOutUMIVsBulk(filtCortRemRUF1, loess1)
regrOut1Lin = regrOutUMIVsBulk(filtCortRemRUF1, lm1)

plot(filtCortRemRUF1$remUMIFrac, filtCortRemRUF1$LogUMIDivBulk)
abline(lm1,col="red")

plot(regrOut1Loess$remUMIFrac, regrOut1Loess$LogUMIDivBulk)
lm2 <- lm(form1, regrOut1Loess)
abline(lm2,col="green") # should be reasonably flat

plot(regrOut1Lin$remUMIFrac, regrOut1Lin$LogUMIDivBulk)
lm3 <- lm(form1, regrOut1Lin)
abline(lm3,col="blue") # should be completely flat

#All looks fine!

