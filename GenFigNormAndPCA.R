#Generates Fig 1 and 2 in the paper.
#Run GenerateMixedDataMatrices before running this file!
library("ggplot2")
library("Seurat")

dataFolder = "C:/Work/R/RNASeqCTProfileEval/"
source(paste0(dataFolder, "FigureHelpFunc.R"))


#############################
# load data
#############################
library(dplyr)
library(tibble)

tpmScAndBulk = readRDS(paste0(dataFolder, "data/tpmScAndBulk.RDS"))
tmmScAndBulk = readRDS(paste0(dataFolder, "data/tmmScAndBulk.RDS"))
quantileScAndBulk = readRDS(paste0(dataFolder, "data/quantileScAndBulk.RDS"))
bcScAndBulk = readRDS(paste0(dataFolder, "data/bcScAndBulk.RDS"))
bcScAndBulkNoCT = readRDS(paste0(dataFolder, "data/bcScAndBulkNoCT.RDS"))


cellTypes = readRDS(paste0(dataFolder, "data/cellTypes.RDS"))
labs = readRDS(paste0(dataFolder, "data/labs.RDS"))
subCellTypes = readRDS(paste0(dataFolder, "data/subCellTypes.RDS"))
tissues = readRDS(paste0(dataFolder, "data/tissues.RDS"))



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

print(paste0("Genes fig 1: ", dim(a)[1]))

library(dplyr)
library(tidyr)

#labvals = c(3,2,4,5,1,6,7,8,9,10)
labvals = c(1,2,3,4,5,6,7,8,9,10)
labnames = c("Bulk 1", "Bulk 2", "Bulk 3", "Bulk 4", "Bulk 5", "SC HCA CB", "SC LC", "SC PBMC68k", "SC Mixed 10x", "SC Melanoma")

#plotSampleOrder = c(2:7,1,8:12,25,13:24,26:64,65:74,75:90,91:98,104:105,101:103,99:100)
plotSampleOrder = c(65:74,1,8:12,25,2:7,13:24,26:64,75:90,91:98,104:105,101:103,99:100)
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

colors = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99",
           "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a") #created with colorbrewer

p1 = ggplot(data= df , aes(x=plotSampleOrder, y=value, fill=Dataset)) +
  geom_boxplot(outlier.shape = NA, coef = 0) +
  labs(y="Relative Log Expression", x="Samples") +
  scale_fill_manual(values=colors) +
  coord_cartesian(ylim=c(-2, 2)) +
  theme_bw() +
  theme( axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
  geom_line(data = dfline, aes(x=xline, y=yline, fill=NA)) + 
  facet_wrap( ~ nm, nrow=3)


#check that the title is shown on the graph, it sometimes randomly disappears. 

library("ggpubr")

fig1 = p1

#fig1 = annotate_figure(p1,
#                top = text_grob("Relative Log Expression After Normalization", face = "bold", size = 14))

fig1

ggsave(
  paste0(fig_path, "Fig1.png"),
  plot = fig1,
  width = 6, height = 4, dpi = 300)

ggsave(
  paste0(fig_path, "Fig1.tiff"),
  plot = fig1,
  width = 6, height = 4, dpi = 300)


###########################
## Figure 2 - PCA plots
###########################

#log transform the data and run PCA
scale = F

a2 = prcomp(t(log2(tpmScAndBulk + 1)), scale=scale)
b2 = prcomp(t(log2(tmmScAndBulk + 1)), scale=scale)
c2 = prcomp(t(log2(quantileScAndBulk + 1)), scale=scale)
d2 = prcomp(t(log2(bcScAndBulk + 1)), scale=scale)

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
#so, here we switch melanoma and PBMC68k to get them in the same alphabetical order (and same colors) as fig 1:
labvals2 = c(1,2,3,4,5,6,7,8,9,10)
labnames2 = c("Bulk 1", "Bulk 2", "Bulk 3", "Bulk 4", "Bulk 5", "SC HCA CB", "SC LC", "SC Melanoma", "SC Mixed 10x", "SC PBMC68k")
labs2 = labs
labs2[labs == 8] = 10
labs2[labs == 10] = 8



dfTPM <- as.data.frame(a2$x)
dfTPM$labs <- factor(labs2, labvals2, labnames2)
dfTPM$bulkAndCt <- factor((cellTypes - 1) + (labs > 5) * 2, c(0,1,2,3), c("Bulk B Cell", "Bulk T Cell", "Sc B Cell", "Sc T Cell"))
dfTPM$pc1Expl = explVar(a2, 1)
dfTPM$pc2Expl = explVar(a2, 2)

dfTMM <- as.data.frame(b2$x)
dfTMM$labs <- factor(labs2, labvals2, labnames2)
dfTMM$bulkAndCt <- factor((cellTypes - 1) + (labs > 5) * 2, c(0,1,2,3), c("Bulk B Cell", "Bulk T Cell", "Sc B Cell", "Sc T Cell"))
dfTMM$pc1Expl = explVar(b2, 1)
dfTMM$pc2Expl = explVar(b2, 2)

dfQ <- as.data.frame(c2$x)
dfQ$labs <- factor(labs2, labvals2, labnames2)
dfQ$bulkAndCt <- factor((cellTypes - 1) + (labs > 5) * 2, c(0,1,2,3), c("Bulk B Cell", "Bulk T Cell", "Sc B Cell", "Sc T Cell"))
dfQ$pc1Expl = explVar(c2, 1)
dfQ$pc2Expl = explVar(c2, 2)

dfBC <- as.data.frame(d2$x)
dfBC$labs <- factor(labs2, labvals2, labnames2)
dfBC$bulkAndCt <- factor((cellTypes - 1) + (labs > 5) * 2, c(0,1,2,3), c("Bulk B Cell", "Bulk T Cell", "Sc B Cell", "Sc T Cell"))
dfBC$pc1Expl = explVar(d2, 1)
dfBC$pc2Expl = explVar(d2, 2)


#Combine all methods:
df = rbind(dfTPM,dfTMM, dfQ, dfBC)
len = dim(dfTPM)[1]
nm = rep(c(1,2,3,4), each=len)
nm = factor(nm, c(1,2,3,4), c("TPM/CPM", "TMM", "Quantile", "TMM + ComBat"))
df = cbind(df,nm)

colors = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99",
           "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a") #created with colorbrewer


p<-ggplot(df,aes(x=PC1,y=PC2,color=labs, shape=bulkAndCt))
p<-p + geom_point()
p<-p + scale_color_manual(values=colors)
p$labels$shape = "Sample prop."
p$labels$colour = "Dataset"
p<-p + coord_cartesian(xlim=c(-0.09, 1), ylim=c(-0.09, 1))#create room for PC label

p<-p + theme_bw() 
p<-p + theme( axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) 

p<-p + geom_rect(data=dat2, mapping=aes(xmin=-0.12, xmax=-0.04, ymin=0, ymax=explPC2), color=NA, fill="blue", alpha=0.2, inherit.aes = FALSE)
p<-p + geom_rect(data=dat2, mapping=aes(xmin=-0.12, xmax=-0.04, ymin=explPC2, ymax=1), color=NA, alpha=0.0, inherit.aes = FALSE)
p<-p + geom_rect(data=dat2, mapping=aes(xmin=-0.12, xmax=-0.04, ymin=0, ymax=1), color="blue", fill=NA, inherit.aes = FALSE)

p<-p + geom_rect(data=dat2, mapping=aes(ymin=-0.124, ymax=-0.04, xmin=0, xmax=explPC1), color=NA, fill="blue", alpha=0.2, inherit.aes = FALSE)
p<-p + geom_rect(data=dat2, mapping=aes(ymin=-0.124, ymax=-0.04, xmin=explPC1, xmax=1), color=NA, alpha=0.0, inherit.aes = FALSE)
p<-p + geom_rect(data=dat2, mapping=aes(ymin=-0.124, ymax=-0.04, xmin=0, xmax=1), color="blue", fill=NA, inherit.aes = FALSE)

p<-p + geom_text(data=dat2, mapping=aes(y=-0.075, x=0.5, label=textPC1), size=3.2, color="black", inherit.aes = FALSE)
p<-p + geom_text(data=dat2, mapping=aes(x=-0.085, y=0.5, label=textPC2), size=3.2, color="black", angle = 90, inherit.aes = FALSE)
p<-p + facet_wrap( ~ nm, nrow=2) 

fig2 = p

#fig2 = annotate_figure(p,
#                top = text_grob("Visualization of Batch Effects", face = "bold", size = 14))
fig2
ggsave(
  paste0(fig_path, "Fig2.png"),
  plot = fig2,
  width = 6, height = 5, dpi = 300)

ggsave(
  paste0(fig_path, "Fig2.tiff"),
  plot = fig2,
  width = 6, height = 5, dpi = 300)

## Figure S1 - PCA plot for combat without cell type in the design matrix
###########################

#log transform the data and run PCA
scale = F

e2 = prcomp(t(log2(bcScAndBulkNoCT + 1)), scale=scale)

explPC1S1=explVar(e2, 1)
explPC2S2=explVar(e2, 2)
textPC1S1=paste0(round(100*explPC1S1,2), "% - Cell type")
textPC2S2=paste0(round(100*explPC2S1,2), "% - Unknown")
datS1 = data.frame(explPC1S1,explPC2S2)

summary(e2) #check that the PC1 explained variance fraction was correctly calculated -> OK!


#scale all PCs to be between 0 and 1
scalePCs <- function(d) {
  for (i in 1:(dim(d)[2])) {
    minval = min(d[,i])
    maxval = max(d[,i])
    d[,i] = (d[,i] - minval) / (maxval - minval)
  }
  return (d)
}
e2$x = scalePCs(e2$x)
#flip PC2 upside down for c so the graphs are more similar. Same for PC1 for b
#c2$x[,2] = -c2$x[,2] + 1
#b2$x[,1] = -b2$x[,1] + 1


#plot(a2$x[,1], a2$x[,2])
#so, here we switch melanoma and PBMC68k to get them in the same alphabetical order (and same colors) as fig 1:
labvals2 = c(1,2,3,4,5,6,7,8,9,10)
labnames2 = c("Bulk 1", "Bulk 2", "Bulk 3", "Bulk 4", "Bulk 5", "SC HCA CB", "SC LC", "SC Melanoma", "SC Mixed 10x", "SC PBMC68k")
labs2 = labs
labs2[labs == 8] = 10
labs2[labs == 10] = 8

dfBCNoCT <- as.data.frame(e2$x)
dfBCNoCT$labs <- factor(labs2, labvals2, labnames2)
dfBCNoCT$bulkAndCt <- factor((cellTypes - 1) + (labs > 5) * 2, c(0,1,2,3), c("Bulk B Cell", "Bulk T Cell", "Sc B Cell", "Sc T Cell"))
dfBCNoCT$pc1Expl = explVar(d2, 1)
dfBCNoCT$pc2Expl = explVar(d2, 2)


#Combine all methods:
#df = rbind(dfTPM,dfTMM, dfQ, dfBC)
#len = dim(dfTPM)[1]
#nm = rep(c(1,2,3,4), each=len)
#nm = factor(nm, c(1,2,3,4), c("TPM/CPM", "TMM", "Quantile", "TMM + ComBat"))
#df = cbind(df,nm)

colors = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99",
           "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a") #created with colorbrewer


p<-ggplot(dfBCNoCT,aes(x=PC1,y=PC2,color=labs, shape=bulkAndCt))
p<-p + geom_point()
p<-p + scale_color_manual(values=colors)
p$labels$shape = "Sample prop."
p$labels$colour = "Dataset"
p<-p + coord_cartesian(xlim=c(-0.09, 1), ylim=c(-0.09, 1))#create room for PC label

p<-p + theme_bw() 
p<-p + theme( axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) 

p<-p + geom_rect(data=datS1, mapping=aes(xmin=-0.12, xmax=-0.04, ymin=0, ymax=explPC2S1), color=NA, fill="blue", alpha=0.2, inherit.aes = FALSE)
p<-p + geom_rect(data=datS1, mapping=aes(xmin=-0.12, xmax=-0.04, ymin=explPC2S1, ymax=1), color=NA, alpha=0.0, inherit.aes = FALSE)
p<-p + geom_rect(data=datS1, mapping=aes(xmin=-0.12, xmax=-0.04, ymin=0, ymax=1), color="blue", fill=NA, inherit.aes = FALSE)

p<-p + geom_rect(data=datS1, mapping=aes(ymin=-0.124, ymax=-0.04, xmin=0, xmax=explPC1S1), color=NA, fill="blue", alpha=0.2, inherit.aes = FALSE)
p<-p + geom_rect(data=datS1, mapping=aes(ymin=-0.124, ymax=-0.04, xmin=explPC1S1, xmax=1), color=NA, alpha=0.0, inherit.aes = FALSE)
p<-p + geom_rect(data=datS1, mapping=aes(ymin=-0.124, ymax=-0.04, xmin=0, xmax=1), color="blue", fill=NA, inherit.aes = FALSE)

p<-p + geom_text(data=datS1, mapping=aes(y=-0.075, x=0.5, label=textPC1S1), size=3.2, color="black", inherit.aes = FALSE)
p<-p + geom_text(data=datS1, mapping=aes(x=-0.085, y=0.5, label=textPC2S1), size=3.2, color="black", angle = 90, inherit.aes = FALSE)
#p<-p + facet_wrap( ~ nm, nrow=2) 

figS1 = p

#fig2 = annotate_figure(p,
#                top = text_grob("Visualization of Batch Effects", face = "bold", size = 14))
figS1
ggsave(
  paste0(fig_path, "FigS1.png"),
  plot = figS1,
  width = 6, height = 5, dpi = 300)




