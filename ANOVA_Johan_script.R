library(xlsx)
library(reshape2)
library(relaimpo)
library(BBmisc)
library(wordspace)
#----------------------- Read in the gene expression data (TMM normalized)
B = read.table(file="C:/Allt/Johan G/bcellProfilesTMMNormalized.txt", skip=1, row.names=1)
head(B)
dim(B)


#------------------ Mean centering
B = log2(B + 0.05)
B = t(B)
is.numeric(B)
B = scale(B, center = T)
B = t(B)
head(B)
sum(B[1, ])
nms = row.names(B)
B = cbind(nms, B)
head(B)

#------------------------ Make into long format
B = melt(B)
head(B)
B = B[order(B[ ,1]), ]
head(B)
B = B[ ,c(1,3)]
head(B)
colnames(B) = c("Gene", "Gene_exp")
head(B)
B[1:35, ]

#----------------------- Read in sample info
meta = read.xlsx(file="C:/Allt/Johan G/SampleCellInfo.xlsx", sheetIndex =1 , startRow=55, row.names=1)
meta = meta[c(1,2,6), ]
meta
names(meta) = NULL
meta
dim(meta)
meta = t(meta)
meta

#----------------------- Merge GE data with sample info
for(i in 1:14) {
  meta = rbind(meta, meta)
}
meta
dim(meta)
dim(B)
meta = meta[1:573370, ]
dim(meta)
B = cbind(B, meta)
head(B)

#-------------------------- REGRESSION 
head(B)
fit = lm(Gene_exp ~ Lab + Tissue + Celltype, data = B)
summary(fit)

lmg = calc.relimp(fit, type = c("lmg", "last", "first", "pratt"), rank = TRUE, diff = TRUE, rela = TRUE)$lmg
last = calc.relimp(fit, type = c("lmg", "last", "first", "pratt"), rank = TRUE, diff = TRUE, rela = TRUE)$last
first = calc.relimp(fit, type = c("lmg", "last", "first", "pratt"), rank = TRUE, diff = TRUE, rela = TRUE)$first
pratt = calc.relimp(fit, type = c("lmg", "last", "first", "pratt"), rank = TRUE, diff = TRUE, rela = TRUE)$pratt
betasq = calc.relimp(fit, type = c("betasq"), rank = TRUE, diff = TRUE, rela = TRUE)$betasq
car = calc.relimp(fit, type = c("car"), rank = TRUE, diff = TRUE, rela = TRUE)$car
genizi = calc.relimp(fit, type = c("genizi"), rank = TRUE, diff = TRUE, rela = TRUE)$genizi

mat = cbind(lmg, last, first, pratt, betasq, car, genizi)
means = c()
sd = c()
for(i in 1:dim(mat)[1]) {
  means = c(means, mean(mat[i,]) )
  sd = c(sd, sd(mat[i,]) )
}
means = as.data.frame(means)
means
dim(means)
row.names(means) = row.names(mat)
# Make better row names
row.names(means) = c("Lab", "Tissue", "Cell type")
means = means[order(-means$means), ,drop = FALSE]
means = t(means)
means
sd

#png(file="C:/Allt/Martin/CIII kinetics paper/Relaimpo_CIII_SR_FCR.png", res=800, height=8, width=5.5, unit = "in")
par(mar=c(8,4,6,6))
barplot(means*100, las=2, ylim=c(0,100), main="Total model adjusted R squared = 0.8 %", 
        ylab="Mean explained variation")
x = 0.68
arrows(x, means[1]*100-sd[1]*100, x, means[1]*100+sd[1]*100, length=0.05, angle=90, code=3, lwd=1)
x = 1.9
arrows(x, means[2]*100-sd[2]*100, x, means[2]*100+sd[2]*100, length=0.05, angle=90, code=3, lwd=1)

#dev.off()



#----------------------------------------- ANOVA
fitaov = aov(Gene_exp ~ Gene + Lab + Tissue + Celltype, data = B[1:35000, ])
ano = anova(fitaov)

summary(ano)




#----------------------------------------------------------------------- Repeated measures ANOVA
library(lme4)
library(lmerTest)
B = read.table(file="C:/Allt/Johan G/bcellProfilesTMMNormalized.txt", skip=1)
head(B)
dim(B)
ngenes = dim(B)[1]
#------------------------ Make into long format
B = melt(B)
head(B)
B = B[order(B[ ,1]), ]
head(B)
B = B[ ,c(1,3)]
head(B)
colnames(B) = c("Gene", "Gene_exp")
head(B)
#----------------------- Read in sample info
meta = read.xlsx(file="C:/Allt/Johan G/SampleCellInfo.xlsx", sheetIndex =1 , startRow=66, row.names=1)
meta = meta[c(1,2,6), ]
meta
names(meta) = NULL
meta
dim(meta)
meta = t(meta)
meta
#----------------------- Merge GE data with sample info
for(i in 1:14) {
  meta = rbind(meta, meta)
}
meta
dim(meta)
dim(B)
meta = meta[1:573370, ]
dim(meta)
B = cbind(B, meta)
head(B)
B[1:36, ]
head(B[-1, ])
head(B)


#---------------------------------- Remove lab 2, 3 and 4 (i.e. row  9-14, (9-14)+35 etc)
a = c()
for (i in 1:ngenes) {
  a = c(a, seq(9,14, by=1)+((i-1)*35))
}
a
length(a)

B = B[-a, ]
dim(B)
B[1:30, ]

#------------------------ filter away lowly expressed genes

d = density(B[ , "Gene_exp"])
plot(d, xlim=c(0,10000))

#---------- log the data
B[ , "Gene_exp"] = log2(B[ , "Gene_exp"] + 0.05)
head(B)
d = density(B[ , "Gene_exp"])
plot(d, xlim=c(-6,2))

B_filt = c()
for (i in 1:ngenes) {
  if( mean(B[(1:35)+((i-1)*35), "Gene_exp"]) > -2 ) {
    B_filt = rbind(B_filt, B[(1:35)+((i-1)*35), ] )
  }
}
dim(B_filt)
head(B_filt)


#-------------------------------- Now run lmer
rmaModel = lmer(Gene_exp ~ Lab + Tissue + Celltype + (1|Gene), data = B)
anova(rmaModel)
coefs <- data.frame(coef(summary(rmaModel)))
# use normal distribution to approximate p-value
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs

#----------------------------------------- ANOVA
psLab = c()
psTissue = c()
psCelltype = c()
for (i in 1:ngenes) {
    if ( mean(B[(1:35)+((i-1)*35), "Gene_exp"]) )
    fitaov = aov(Gene_exp ~  Lab + Tissue + Celltype, data = B[(1:35)+((i-1)*35), ])
    psLab = c(psLab, summary(fitaov)[[1]][["Pr(>F)"]][1] )
    psTissue = c(psTissue, summary(fitaov)[[1]][["Pr(>F)"]][2] )
    psCelltype = c(psCelltype, summary(fitaov)[[1]][["Pr(>F)"]][3] )
  
}
median(psLab, na.rm=T)
median(psTissue, na.rm=T)
median(psCelltype, na.rm=T)

boxplot(psLab)
boxplot(psTissue)
boxplot(psCelltype)

fitaov = aov(Gene_exp ~  Gene + Lab + Tissue + Celltype, data = B[1:100000, ])
summary(fitaov)

#------------------------------------- MY OWN ANALYSIS (HEMSNICKRAD)
B = read.table(file="C:/Allt/Johan G/bcellProfilesTMMNormalized.txt", skip=1, row.names=1)
head(B)
dim(B)

#------------------ Mean centering
B = log2(B + 0.05)
B = t(B)
is.numeric(B)
B = scale(B, center = T)
B = t(B)
head(B)
sum(B[1, ])
nms = row.names(B)
#B = cbind(nms, B)
head(B)

#------------------------ Make into long format
B = melt(B)
head(B)
B = B[order(B[ ,1]), ]
head(B)
#B = B[-1, ]
B = B[ ,c(1,3)]
head(B)
colnames(B) = c("Gene", "Gene_exp")
head(B)
B[1:37, ]

#----------------------- Read in sample info
meta = read.xlsx(file="C:/Allt/Johan G/SampleCellInfo.xlsx", sheetIndex =1 , startRow=55, row.names=1)
meta = meta[c(1,2,6), ]
meta
names(meta) = NULL
meta
dim(meta)
meta = t(meta)

#----------------------- Merge GE data with sample info
for(i in 1:14) {
  meta = rbind(meta, meta)
}
meta
dim(meta)
dim(B)
meta = meta[1:573370, ]
dim(meta)
B = cbind(B, meta)
head(B)
B[1:35,]

#--------------- Loop over data and calculate means for every gene
means_tot = c()
means_Lab1 = c()
means_Lab2 = c()
means_Lab3 = c()
means_Lab4 = c()
means_Lab5 = c()
means_Lab7 = c()
means_Tissue1 = c()
means_Tissue2 = c()
means_Tissue3 = c()
means_Tissue4 = c()
means_Celltype1 = c()
means_Celltype2 = c()
means_Celltype4 = c()
means_Celltype5 = c()
means_Celltype6 = c()
means_Celltype8 = c()


for ( i in 1:ngenes) {
  means_tot = c(means_tot, mean(B[(1:35)+((i-1)*35), "Gene_exp"]) )
  
  means_Lab1 = c(means_Lab1, mean(B[(1:8)+((i-1)*35), "Gene_exp"]) )
  means_Lab2 = c(means_Lab2, mean(B[(14)+((i-1)*35), "Gene_exp"]) )
  means_Lab3 = c(means_Lab3, mean(B[(12:13)+((i-1)*35), "Gene_exp"]) )
  means_Lab4 = c(means_Lab4, mean(B[(9:11)+((i-1)*35), "Gene_exp"]) )
  means_Lab5 = c(means_Lab5, mean(B[(15:30)+((i-1)*35), "Gene_exp"]) )
  means_Lab7 = c(means_Lab7, mean(B[(31:35)+((i-1)*35), "Gene_exp"]) )

  means_Tissue1 = c(means_Tissue1, mean(B[(which(B[1:35, "Tissue"] == "A"))+((i-1)*35), "Gene_exp"]) )
  means_Tissue2 = c(means_Tissue2, mean(B[(which(B[1:35, "Tissue"] == "B"))+((i-1)*35), "Gene_exp"]) )
  means_Tissue3 = c(means_Tissue3, mean(B[(which(B[1:35, "Tissue"] == "C"))+((i-1)*35), "Gene_exp"]) )
  means_Tissue4 = c(means_Tissue4, mean(B[(which(B[1:35, "Tissue"] == "D"))+((i-1)*35), "Gene_exp"]) )
  
  means_Celltype1 = c(means_Celltype1, mean(B[(which(B[1:35, "Celltype"] == "A"))+((i-1)*35), "Gene_exp"]) )
  means_Celltype2 = c(means_Celltype2, mean(B[(which(B[1:35, "Celltype"] == "B"))+((i-1)*35), "Gene_exp"]) )
  means_Celltype4 = c(means_Celltype4, mean(B[(which(B[1:35, "Celltype"] == "D"))+((i-1)*35), "Gene_exp"]) )
  means_Celltype5 = c(means_Celltype5, mean(B[(which(B[1:35, "Celltype"] == "E"))+((i-1)*35), "Gene_exp"]) )
  means_Celltype6 = c(means_Celltype6, mean(B[(which(B[1:35, "Celltype"] == "F"))+((i-1)*35), "Gene_exp"]) )
  means_Celltype8 = c(means_Celltype8, mean(B[(which(B[1:35, "Celltype"] == "H"))+((i-1)*35), "Gene_exp"]) )

}
mean(means_tot, na.rm=T)

par(mfrow=c(1,3))
nLab1 = 8
nLab2 = 1
nLab3 = 
boxplot(means_tot, means_Lab1, means_Lab2, means_Lab3, means_Lab4, means_Lab5, means_Lab7, main = "Lab")
boxplot(means_tot, means_Tissue1, means_Tissue2, means_Tissue3, means_Tissue4, main = "Tissue")
boxplot(means_tot, means_Celltype1, means_Celltype2, means_Celltype4, means_Celltype5, means_Celltype6, means_Celltype8, main = "Cell type")






