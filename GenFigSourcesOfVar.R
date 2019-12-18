#Run Generate MixedDataMatrices before running this file!
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

tpmScAndBulkNonFilt = readRDS(paste0(dataFolder, "data/tpmScAndBulkNonFilt.RDS"))
tmmScAndBulkNonFilt = readRDS(paste0(dataFolder, "data/tmmScAndBulkNonFilt.RDS"))
bcScAndBulk = readRDS(paste0(dataFolder, "data/quantileScAndBulkNonFilt.RDS"))

cellTypes = readRDS(paste0(dataFolder, "data/cellTypes.RDS"))
labs = readRDS(paste0(dataFolder, "data/labs.RDS"))
subCellTypes = readRDS(paste0(dataFolder, "data/subCellTypes.RDS"))
tissues = readRDS(paste0(dataFolder, "data/tissues.RDS"))
individual = readRDS(paste0(dataFolder, "data/individual.RDS"))
techRepl = readRDS(paste0(dataFolder, "data/techRepl.RDS"))
scOrBulk = readRDS(paste0(dataFolder, "data/scOrBulk.RDS"))


###########################
## Special calculation for individual and technical replicates
###########################
labIds = unique(labs)
#figure out how many within-lab combinations there can max be
maxSamplesPerLab = max(table(labs)) #table gives samples per lab
maxComb = maxSamplesPerLab * (maxSamplesPerLab-1)/2
samps = 1:length(labs)
pseudoCount = 0.05
tecRep = matrix(data=NA, nrow=maxComb, ncol=length(labIds))
sameInd = matrix(data=NA, nrow=maxComb, ncol=length(labIds))
diffInd = matrix(data=NA, nrow=maxComb, ncol=length(labIds))
#so, these are matrices where each column represents a lab
#the matrices are just filled with values from pairs, as many as are found.

#loop through all labs
for (lab in labIds) {
  #keep track of how the next free index in each vector
  ti = 1
  si = 1
  di = 1
  #loop through all combinations
  for (samp1 in 1:(length(labs)-1)) {
    if (labs[samp1] == lab) {
      for (samp2 in (samp1+1):length(labs)) {
        #check that we're in the same lab and in the same tissue, cell type and sub cell type, otherwise just ignore the pair
        if ((labs[samp2] == lab) & 
            (cellTypes[samp1] == cellTypes[samp2]) &
            (subCellTypes[samp1] == subCellTypes[samp2]) &
            (tissues[samp1] == tissues[samp2])
            ) {
          #calculate difference
          diff = sd(log2((tmmScAndBulk[,samp1] + pseudoCount)/(tmmScAndBulk[,samp2] + pseudoCount)))
          #figure out if it is a technical replicate, same individual, or different individuals
          #technical replicate
          if ((!is.na(techRepl[samp1])) & (!is.na(techRepl[samp2])) & ( techRepl[samp1] == techRepl[samp2] ) ) {
            tecRep[ti,lab] = diff
            ti = ti+1
          }
          #same individual, but not technical replicate
          else if ((!is.na(individual[samp1])) & (!is.na(individual[samp2])) & ( individual[samp1] == individual[samp2] ) ) {
            sameInd[si,lab] = diff
            si = si+1
          }
          #different individuals
          else {
            diffInd[di,lab] = diff
            di = di+1
          }
        }
      }
    }
  }
}

#It seems that values for both technical replicates and different individuals are 
#only present for patient 2. Calculate the means of those two groups:
techRepLab3 = tecRep[!is.na(tecRep[,3]),3]
diffLab3 = diffInd[!is.na(diffInd[,3]),3]
techRepfact = mean(techRepLab3)/mean(diffLab3)

sameIndLab4 = sameInd[!is.na(sameInd[,4]),4]
diffIndLab4 = diffInd[!is.na(diffInd[,4]),4]
sameIndfact = mean(sameIndLab4)/mean(diffIndLab4)



library(ggplot2)

boxPlotWithDots <- function(data, boxes, color){
  
  #randomize x:es:
  xes = runif(length(data)) - 0.5
  df <- data.frame(d=data, x=xes, boxes = boxes)

  g1 <- ggplot(df, aes(y=d))+geom_boxplot(outlier.shape = NA) +
    geom_point(alpha=0.5, aes(x=x, color=boxes),size=2) + 
    labs(title="Technical Replicates and Individuals", y="Std(log fold change)", x="")
  
  g1 = g1+theme(panel.grid.major= element_blank(),
                panel.grid.minor= element_blank(),
                panel.background= element_blank(),
                panel.border= element_blank(),
                legend.position="none",
                axis.text.x=element_blank(), 
                axis.ticks.x=element_blank(), 
                axis.ticks.y=element_blank())
  #g1= g1+scale_shape_manual(values=c(19,1))+scale_fill_discrete(guide=FALSE)
#  print(g1)
}

joinedData = c(techRepLab3, diffLab3, sameIndLab4, diffIndLab4)
boxes = c(rep(1,length(techRepLab3)), rep(2,length(diffLab3)), rep(3,length(sameIndLab4)), rep(4,length(diffIndLab4)))
boxFac = factor(boxes, c(1,2,3,4), c("Techn. repl. bulk 3", "Diff. indiv. bulk 3", "Same indiv. bulk 4", "Diff. indiv. bulk 4"))


p = boxPlotWithDots(joinedData, boxFac, 1)
p<-p + facet_wrap( ~ boxes, nrow=1) 
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
