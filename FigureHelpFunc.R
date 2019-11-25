
#############################
# Help functions
#############################
MakeTPM <- function(ds) {
  cs = colSums(ds)
  for(i in 1:length(cs)) {
    ds[,i] = ds[,i]/cs[i] * 10^6
  }
  return (ds)
}

library("edgeR")

TMMNorm <- function(ds) {
  #using TMM from edgeR:
  normFactors <- calcNormFactors(ds,NULL,"TMM")
  
  libSizes = colSums(ds)
  effectiveLibSizes = libSizes * normFactors
  
  
  #I need to transpose the data matrix back and forth to get the
  #row wise division to work...
  ds = t(t(ds) * (10^6/ effectiveLibSizes))
  
  #There are some issues with roundoff or similar making the mean of the sum of all 
  #genes not to be exactly 10^6. Fix this:
  sumPerSamp = colSums(ds)
  meanSampSum = mean(sumPerSamp)
  
  ds = ds * 10^6/meanSampSum
  
  
  return (ds)
}


#to be skipped later
GetScBias <- function(expr1, expr2, lb=-3, ub=12) {
  divvar = log2((expr1 + 0.05) / (expr2 + 0.05))
  #meanExprVar = log2((expr1 + expr2)/2 + 0.05)
  meanExprVar = log2(sqrt(expr1*expr2) + 0.05)
  #meanExprVar = log2(expr2 + 0.05)
  
  step = 0.5
  aimedxes = seq(lb, ub, by=step)
  uppers = aimedxes + step/2
  lowers = aimedxes - step/2
  
  #allocate of the correct length
  x = aimedxes
  y = aimedxes
  
  for (i in 1:length(aimedxes)) {
    ind = (meanExprVar >= lowers[i]) & (meanExprVar <= uppers[i])
    x[i] = mean(meanExprVar[ind])
    y[i] = mean(divvar[ind])
  }
  
  return (cbind(x, y))
}


corrUMIVsBulk <- function(ds) {
  cor(ds$logUMITMM, ds$logBulkTMM)
}

#will regress out in log space and update the fields "logUMITMM" and "LogUMIDivBulk"
#may produce some warnings about missing values that can be ignored
regrOutUMIVsBulk <- function(ds, fit) {
  ds2 = ds
  pred = predict(fit,ds2)
  fitNAFilter = !is.na(pred);
  #Now, regress out
  ds2$logUMITMM[fitNAFilter] = ds2$logUMITMM[fitNAFilter] - pred[fitNAFilter] + mean(pred[fitNAFilter])
  #We should not restore the genes that has become really lowly expressed after this many of them
  #are lowly expressed in bulk
  #ds2$logUMITMM[ds2$logUMITMM < log2(1.05)] = log2(1.05) #the data is filtered on 1 TPM ~ 1 in TMM, so no values should be below that
  
  #also update the "LogUMIDivBulk"
  ds2$LogUMIDivBulk = ds2$logUMITMM - ds2$logBulkTMM;
  
  return (ds2)
}

genScToBulkCovGraphs <- function(ds, formulaUMI, formulaLFC, covIndex, covName, filter = NA) {
  if (!is.na(filter[1])) {
    ds = ds[filter,]
    print("filtering")
  }
  ind = sort(ds[,covIndex], index.return=T, na.last = T)
  dsSort = ds[ind$ix,];
  naFilt = !is.na(dsSort[,covIndex])
  numGenes = dim(dsSort)[1]
  dsPlot = data.frame(dsSort[,covIndex], dsSort$logUMITMM)
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
  #lines(dsSort[,covIndex], predict(lm2), col = "blue")
  
  #Plot log fold change between UMI and bulk vs covariate 
  dsPlot = data.frame(dsSort[,covIndex], dsSort$LogUMIDivBulk)
  colnames(dsPlot) = c("x","y")
  
  dsLoess = data.frame(dsSort[naFilt,covIndex], predict(loess_fit2, dsSort[naFilt,]))
  colnames(dsLoess) = c("x","y")
  dsLin = data.frame(dsSort[naFilt,covIndex], predict(lm2, dsSort[naFilt,]))
  colnames(dsLin) = c("x","y")
  
  p2 = ggplot(dsPlot,aes(x=x,y=y))
  p2 = p2 + geom_point(alpha=0.3, shape=1)
  p2 = p2 + geom_line(data = dsLoess, colour="#FF0000", alpha=1, size=1.4)
  p2 = p2 + geom_line(data = dsLin, colour="#0000FF", alpha=1, size=1.4)
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
  print(ruvsb)
  ruvsbl = corrUMIVsBulk(dsRegrLin)
  print(ruvsbl)
  #plot(dsRegr[,covIndex], dsRegr$LogUMIDivBulk) # for test only

  res = list(p1, p2, uvsb, ruvsb, ruvsbl)
  
  return (res)
}

#assumes the following structure: Bulk1 Bulk2 BulkCounts1 BulkCounts2 UMI1 UMI2 count1 count2 gcFullLength    gcTail   tx_len
extractSample <- function(mergedData, index) {
  addN = index - 1
  dat = mergedData[,c(1+addN, 5+addN, 9, 10, 11)]
  #calculate removed counts' fraction
  dat = cbind(dat, (mergedData[,7+addN]- mergedData[,5+addN]) / mergedData[,7+addN])
  #calculate removed counts' fraction for the other cortex:
  invAddN = 1-addN
  dat = cbind(dat, (mergedData[,7+invAddN]- mergedData[,5+invAddN]) / mergedData[,7+invAddN])
  
  #add TPM, and TMM-normalized, log transformed data
  d2 = dat[,c(1,2)]
  d2 = MakeTPM(d2);
  #scale according to counts before TMM
  csCounts = sum(mergedData[,3+addN])
  pseudoCounts  = cbind(dat[,1]*(2*csCounts/10^6), dat[,2])
  tmmNorm = TMMNorm(pseudoCounts)
  d3 = log2(tmmNorm + 0.05)
  d3 = cbind(d3, log2((tmmNorm[,2] + 0.05)/(tmmNorm[,1] + 0.05)))
  dat = cbind(dat, d2,d3)
  colnames(dat) = c("bulk", "UMI", "gcFullLength", "gcTail", "geneLength", "remUMIFrac", "remUMIFracOtherSample", "bulkTPM", "UMITPM", "logBulkTMM", "logUMITMM", "LogUMIDivBulk")
  return (dat)
}
