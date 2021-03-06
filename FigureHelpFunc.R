#############################
# Help functions
#############################

#Normalizes to TPM/CPM
MakeTPM <- function(ds) {
  cs = colSums(ds)
  for(i in 1:length(cs)) {
    ds[,i] = ds[,i]/cs[i] * 10^6
  }
  return (ds)
}

#Does TMM normalization using edgeR
TMMNorm <- function(ds) {
  #using TMM from edgeR:
  normFactors <- edgeR::calcNormFactors(ds,NULL,"TMM")
  
  libSizes = colSums(ds)
  effectiveLibSizes = libSizes * normFactors
  
  #I need to transpose the data matrix back and forth to get the
  #row wise division to work...
  ds = t(t(ds) * (10^6/ effectiveLibSizes))
  
  #There are some issues with roundoff or similar making the mean of the sum of all 
  #genes not to be exactly 10^6. Fix this:
  sumPerSamp = colSums(ds)
  meanSampSum = mean(sumPerSamp)
  
  ds = ds * 10^6 / meanSampSum
  
  return (ds)
}


#Help function used from GenFigScVsBulk.
#Creates plot A and B in fig 5.
plotCorr <- function(expr10x, exprBulk) {
  ds = data.frame(exprBulk, expr10x)
  ind = sort(ds[,1], index.return=T, na.last = T)
  dsSort = ds[ind$ix,];
  colnames(dsSort) = c("x","y")
  #loess_fit <- loess(y ~x, dsSort, span = 0.3)
  
  #dsLoess = data.frame(dsSort$x, predict(loess_fit, dsSort)) 
  #colnames(dsLoess) = c("x","y")
  dsProp = data.frame(x=c(-5,15), y=c(-5,15))
  
  #Plot log expression in UMI data vs covariate 
  p1 = ggplot(dsSort,aes(x=x,y=y))
  #p1 = p1 + geom_hex(bins = 70)#geom_point(alpha=0.3, shape=1) 
  p1 = p1 + stat_binhex(aes(fill=log(..count..)), bins = 70)
  #p1 = p1 + scale_fill_continuous(type = "viridis")
  #p1 = p1 + geom_line(data = dsProp, colour="#FF0000", size=1.4, alpha=1)
  p1 = p1 + geom_line(data = dsProp, colour="#FF8800", size=1.4, alpha=1) #good
  #p1 = p1 + geom_line(data = dsProp, colour="#66CC00", size=1.4, alpha=1)
  #p1 = p1 + geom_line(data = dsLoess, colour="#FF0000", size=1.4, alpha=1) # the loess fit doesn't look that good, skip it
  p1 = p1 + labs(y="10X, log2(pseudo-CPM)", x="Bulk, log2(pseudo-TPM)")
  #p1<-p1 + labs(title="Visualization of Batch Effects")
  p1 = p1 + coord_cartesian(xlim=c(-5, 15), ylim=c(-5, 15))#create room for PC label
  #p<-p + theme( axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) 
  p1 = p1 + theme_bw() + theme(legend.position="bottom")
  print(p1)
  
  return (p1)  
}

