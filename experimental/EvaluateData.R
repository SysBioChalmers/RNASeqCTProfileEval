
#First read TPM and count matrix
folder = "C:/Work/R/RNASeqCTProfileEval"
bulk_tpm = read.table(file=paste0(folder, "/tpmMatrix.txt"), header=T, sep="\t")
bulk_tmm = read.table(file=paste0(folder, "/tmmMatrix.txt"), header=T, sep="\t")
bulk_counts = read.table(file=paste0(folder, "/countsMatrix.txt"), header=T, sep="\t")

totCounts = colSums(bulk_counts);

#read the design matrix
#install.packages("xlsx")
library("xlsx")
desm <- read.xlsx(paste0(folder, "/DesignMatrix.xlsx"), sheetName = "DesignMatrix")
cellTypes = as.numeric(desm[9, 2:75])
labs = as.numeric(desm[2, 2:75])
subCellTypes = as.numeric(desm[7, 2:75])
tissues = as.numeric(desm[3, 2:75])


#calculates the variable index for a combination between two samples
#a and b are the variable values from the two samples (such as lab)
#maxval how many labs there are (the max value of the lab variable, which should be the same).
#0 means it is part of intercept
var.index <- function(a, b, maxval) {
  r = min(a,b);
  c = max(a,b);
  #The row is the smallest value, so first add up the rows above
  above = 0;
  if (r > 1) {
    for (i in 1:(r-1)) {
      above = above + (maxval - (i-1));
    }
  }
  #return the result
  above + c - r + 1 -1 #+ 1 to return from 0 indexing, -1 to make the first part of the intercept
}

#gets the sample pair index from two sample indices
get.ind <- function(i1, i2, numSamp) {
  if (i1 == i2) {
    return(0);
  }
  r = min(i1,i2);
  c = max(i1,i2);
  ind = 0;
  #add all rows above i1
  for (i in 1:(r-1)) {
    ind = ind + numSamp - i;
  }
  #add the appropriate number for i2
  return(ind + c-r);
}


#test
var.index(2,3,4) # should be 6
var.index(3,2,4) # should be 6
var.index(6,4,6) # should be 18

get.ind(3,4,6) # should be 10
get.ind(1,1,6) # should be 0
get.ind(4,6,6) # should be 14



#create the regression combinations
#Let's use lab, celltype and tissue as parameters in the first regression
#Let's code it like this: highest number * (number of first-1) + number of second - 1;
#always select the smaller of them as first

#generate the design matrix
numSamp = length(cellTypes);
#numSamp = 10;#keep it short for now

numComb = numSamp*(numSamp-1)/2

numLabs = length(unique(labs))
numCT = length(unique(cellTypes))
numTis = length(unique(tissues))

numLabVars = numLabs*(numLabs+1)/2-1;
numCTVars = numCT*(numCT+1)/2-1;
numTisVars = numTis*(numTis+1)/2-1;

numVars = numLabVars + numCTVars + numTisVars;

dm = as.data.frame(matrix(0, numComb, numVars))
counts1 = rep(0,numComb)
counts2 = rep(0,numComb)


labInd0 = 0;
ctInd0 = labInd0 + numLabVars;
tisInd0 = ctInd0 + numCTVars;


#set the names of the vars
for (i in 1:numLabs) {
  #this will run some things twice, but it matters not
  for (j in 1:numLabs) {
    ind = var.index(i,j,numLabs);
    if (ind != 0) {
      low = min(i,j)
      high = max(i,j)
      colnames(dm)[ind] = paste0("L_",low,"vs",high);
    }
  }
}
for (i in 1:numCT) {
  #this will run some things twice, but it matters not
  for (j in 1:numCT) {
    ind = var.index(i,j,numCT);
    if (ind != 0) {
      low = min(i,j)
      high = max(i,j)
      colnames(dm)[ind + ctInd0] = paste0("CT_",low,"vs",high);
    }
  }
}
for (i in 1:numTis) {
  #this will run some things twice, but it matters not
  for (j in 1:numTis) {
    ind = var.index(i,j,numTis);
    if (ind != 0) {
      low = min(i,j)
      high = max(i,j)
      colnames(dm)[ind + tisInd0] = paste0("TIS_",low,"vs",high);
    }
  }
}



#select the data to work with and filter appropriately
#datMat = bulk_tpm;
datMat = bulk_tmm;

meanExpr = rowSums(datMat) / numSamp;
filt = meanExpr > 1
datMat = datMat[filt,]
y = rep(0, numComb)

numGenes = dim(datMat)[1]
combMat = matrix(0,numGenes, numComb)
combMeanExpr = matrix(0,numGenes, numComb)
dim(combMat)
meanExprFilt = rowMeans(datMat);

#now build the design matrix and dependent variable
index = 1;
for (i in 1:(numSamp-1)) {
  for (j in (i+1):numSamp) {
#    print(index)
#    print (i)
#    print (j)

    #design matrix
    #lab
    ind = var.index(labs[i], labs[j], numLabs)
    if (ind != 0) {
      ii = labInd0 + ind;
      dm[index,ii] = 1;
    }
    #ct
    ind = var.index(cellTypes[i], cellTypes[j], numCT)
    if (ind != 0) {
      ii = ctInd0 + ind;
      dm[index,ii] = 1;
    }
    #tissue
    ind = var.index(tissues[i], tissues[j], numTis)
    if (ind != 0) {
      ii = tisInd0 + ind;
      dm[index,ii] = 1;
    }

    #The data matrix
    a = log2((datMat[,i] + 1)/(datMat[,j] + 1))
    combMat[,index] = a;
    y[index] = sd(a);
    combMeanExpr[,index] = ((datMat[,i]) + (datMat[,j]))/2;

    #tot counts per sample
    counts1[index] = totCounts[i];
    counts2[index] = totCounts[j];

    index = index+1;
  }
}

#now run the regression
dataset = cbind(y,dm)

res = lm(y~.,dataset)

summary(res)


#investigate the relationship between noise and counts
#plot(counts1,y, log="x")
#plot(counts1 * counts2,y, log="x")
#x = log(1/counts1 + 1/counts2)
#plot(x,y, log="x")
library(ggplot2)

#df = as.data.frame(cbind(x,y))

#ggplot(df, aes(x = x, y = y)) +
#  geom_smooth()
#ggplot(1/counts1 + 1/counts2,y, aes(x = Time, y = I)) +
#  geom_point(aes(color = ID)) +
#  geom_smooth()

#hmm, not a super clear relationship. The number of reads are probably
#well correlated with lab, so we probably do not need to take care of that.
#check reads vs lab correlation:
#plot(labs, totCounts, log="y")
#this is absolutely the case


#check the difference between some samples
#technical replicates vs different labs
a = log2((datMat[,3] + 1)/(datMat[,4] + 1))
b = log2((datMat[,17] + 1)/(datMat[,28] + 1))
d <- density(a) # returns the density data
plot(d) # plots the results

c = c(a,b)

df = as.data.frame(c)
#ggplot(data=df, aes(df$a)) +
#  geom_histogram(binwidth = 0.05) +
#  geom_histogram(aes = b, binwidth = 0.05)

num = length(a);
sel = c(rep(0,num),rep(1,num));

ggplot(data=df, aes(x=c)) +
  geom_histogram(data=subset(df,sel == 0),fill = "red", alpha = 0.2, binwidth = 0.05) +
  geom_histogram(data=subset(df,sel == 1),fill = "blue", alpha = 0.2, binwidth = 0.05)


#ggplot(df, aes(x)) +                    # basic graphical object
#  geom_line(aes(y=y1), colour="red") +  # first layer
#  geom_line(aes(y=y2), colour="green")  # second layer

#dat <- data.frame(xx = c(runif(100,20,50),runif(100,40,80),runif(100,0,30)),yy = rep(letters[1:3],each = 100))

#ggplot(dat,aes(x=xx)) +
#  geom_histogram(data=subset(dat,yy == 'a'),fill = "red", alpha = 0.2) +
#  geom_histogram(data=subset(dat,yy == 'b'),fill = "blue", alpha = 0.2) +
#  geom_histogram(data=subset(dat,yy == 'c'),fill = "green", alpha = 0.2)

#same lab vs different labs
a = log2((datMat[,2] + 1)/(datMat[,5] + 1))
b = log2((datMat[,17] + 1)/(datMat[,28] + 1))
#d <- density(a) # returns the density data
#plot(d) # plots the results

c = c(a,b)

df = as.data.frame(c)

num = length(a);
sel = c(rep(0,num),rep(1,num));

ggplot(data=df, aes(x=c)) +
  geom_histogram(data=subset(df,sel == 0),fill = "red", alpha = 0.2, binwidth = 0.05) +
  geom_histogram(data=subset(df,sel == 1),fill = "blue", alpha = 0.2, binwidth = 0.05)


#same individual vs different individuals, but same lab
a = log2((datMat[,13] + 1)/(datMat[,14] + 1))
b = log2((datMat[,13] + 1)/(datMat[,16] + 1))

c = c(a,b)

df = as.data.frame(c)

num = length(a);
sel = c(rep(0,num),rep(1,num));

ggplot(data=df, aes(x=c)) +
  geom_histogram(data=subset(df,sel == 0),fill = "red", alpha = 0.2, binwidth = 0.05) +
  geom_histogram(data=subset(df,sel == 1),fill = "blue", alpha = 0.2, binwidth = 0.05)

#one vs two factors
a = log2((datMat[,28] + 1)/(datMat[,29] + 1))
b = log2((datMat[,15] + 1)/(datMat[,29] + 1))

c = c(a,b)

df = as.data.frame(c)

num = length(a);
sel = c(rep(0,num),rep(1,num));

ggplot(data=df, aes(x=c)) +
  geom_histogram(data=subset(df,sel == 0),fill = "red", alpha = 0.2, binwidth = 0.05) +
  geom_histogram(data=subset(df,sel == 1),fill = "blue", alpha = 0.2, binwidth = 0.05)

#cell type
a = log2((datMat[,17] + 1)/(datMat[,18] + 1))
b = log2((datMat[,18] + 1)/(datMat[,19] + 1))

c = c(a,b)

df = as.data.frame(c)

num = length(a);
sel = c(rep(0,num),rep(1,num));

ggplot(data=df, aes(x=c)) +
  geom_histogram(data=subset(df,sel == 0),fill = "red", alpha = 0.2, binwidth = 0.05) +
  geom_histogram(data=subset(df,sel == 1),fill = "blue", alpha = 0.2, binwidth = 0.05)



# this is how to write "without intercept": ~.-1

#use variance instead of stddev, residuals vs fitted looks much better then, and it would also
#make sense if the factors are multiplicative
variance = y^2


#investigate the linearity
#select pairs of samples that have lab 4,4 or 4, 5, tissue 1, cell type 1, 1 or 1,2
#so, 4 groups of samples

selL4CT1 = which((labs == 4) & (tissues == 1) & (subCellTypes == 1) & (cellTypes == 1))
selL5CT1 = which((labs == 5) & (tissues == 1) & (subCellTypes == 1) & (cellTypes == 1))
selL4CT2 = which((labs == 4) & (tissues == 1) & (subCellTypes == 1) & (cellTypes == 2))
selL5CT2 = which((labs == 5) & (tissues == 1) & (subCellTypes == 3) & (cellTypes == 2))


selL44CT11 = logical(numComb);
selL45CT11 = logical(numComb);
selL44CT12 = logical(numComb);
selL45CT12 = logical(numComb);
#selL44CT11
for (i in selL4CT1) {
  for (j in selL4CT1) {
    ind = get.ind(i,j,numSamp)
    if (ind != 0) {
      selL44CT11[ind] = TRUE;
    }
  }
}

#selL45CT11
for (i in selL4CT1) {
  for (j in selL5CT1) {
    ind = get.ind(i,j,numSamp)
    if (ind != 0) {
      selL45CT11[ind] = TRUE;
    }
  }
}

#selL44CT12
for (i in selL4CT1) {
  for (j in selL4CT2) {
    ind = get.ind(i,j,numSamp)
    if (ind != 0) {
      selL44CT12[ind] = TRUE;
    }
  }
}

#selL45CT12
for (i in selL4CT1) {
  for (j in selL5CT2) {
    ind = get.ind(i,j,numSamp)
    if (ind != 0) {
      selL45CT12[ind] = TRUE;
    }
  }
}

#transvar = variance^2.3
transvar = variance

vL44CT11 = mean(transvar[selL44CT11])
vL45CT11 = mean(transvar[selL45CT11])
vL44CT12 = mean(transvar[selL44CT12])
vL45CT12 = mean(transvar[selL45CT12])

vL44CT11
vL45CT11
vL44CT12
vL45CT12

dataset3 = cbind(transvar,dm)

res3 = lm(transvar~.,dataset3)
par(mfrow=c(2,2))
plot(res3)


#cellTypes = as.numeric(desm[9, 2:75])
#labs = as.numeric(desm[2, 2:75])
#subCellTypes = as.numeric(desm[7, 2:75])
#tissues

par(mfrow=c(2,2))
plot(transvar[selL44CT11])
plot(transvar[selL45CT11])
plot(transvar[selL44CT12])
plot(transvar[selL45CT12])

par(mfrow=c(1,1))
hist(transvar)
summary(res3)

genePool = 300

#make a loop over the gene expresson range (sliding window) to show this as a function of gene expression range:
L44CT11s = rep(0,numGenes - genePool);
L44CT11sx = L44CT11s;
L45CT11s = L44CT11s;
L45CT11sx = L44CT11s;
L44CT12s = L44CT11s;
L44CT12sx = L44CT11s;
L45CT12s = L44CT11s;
L45CT12sx = L44CT11s;

sortedGeneInd = sort(meanExprFilt, index.return=T)
for (i in seq(1, numGenes - genePool, by=50)) {
#  print(i)
  geneSubsetData = combMat[sortedGeneInd$ix[i:(i+genePool-1)],]
  geneSubsetDatMat = combMeanExpr[sortedGeneInd$ix[i:(i+genePool-1)],]
  tempVar = apply(geneSubsetData,2,var);

  L44CT11s[i] = mean(tempVar[selL44CT11])
  L44CT11sx[i] = mean(rowMeans(geneSubsetDatMat[,selL44CT11]))
  L45CT11s[i] = mean(tempVar[selL45CT11])
  L45CT11sx[i] = mean(rowMeans(geneSubsetDatMat[,selL45CT11]))
  L44CT12s[i] = mean(tempVar[selL44CT12])
  L44CT12sx[i] = mean(rowMeans(geneSubsetDatMat[,selL44CT12]))
  L45CT12s[i] = mean(tempVar[selL45CT12])
  L45CT12sx[i] = mean(rowMeans(geneSubsetDatMat[,selL45CT12]))

}

ggplot() + geom_line(aes(x=L44CT11sx,y=L44CT11s),color='red') +
  geom_line(aes(x=L45CT11sx,y=L45CT11s),color='blue') +
  geom_line(aes(x=L44CT12sx,y=L44CT12s),color='green') +
  geom_line(aes(x=L45CT12sx,y=L45CT12s),color='orange') +
  ylab('var')+xlab('tpm')

#check if the fact that we have the same patient matters and makes the red curve too low, if that creates
#a problem here


df = as.data.frame(cbind(L44CT11sx,L44CT11s,L45CT11s,L45CT11sx,L44CT12s,L44CT12sx,L45CT12s,L45CT12sx))




library("reshape2")

test_data_long <- melt(df, id="date")  # convert to long format

ggplot(data=test_data_long,
       aes(x=date, y=value, colour=variable)) +
  geom_line()


ggplot(df, aes(x = L44CT11sx, y = L44CT11s)) + geom_smooth() +
  geom_plot(aes(x = L45CT11sx, y = L45CT11s)) + geom_smooth()

plot(L44CT11sx,L44CT11s)
plot(L44CT11sx,L44CT11s)
plot(L44CT11sx,L44CT11s)
plot(L44CT11sx,L44CT11s)


dataset2 = cbind(variance,dm)

res2 = lm(variance~.,dataset2)
plot(res2)

summary(res2)




