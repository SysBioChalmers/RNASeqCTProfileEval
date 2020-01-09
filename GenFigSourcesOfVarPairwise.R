#This script generates Fig 3 and Fig S1.

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
## Set up all pairs with only one factor
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
diffTissue = matrix(data=NA, nrow=maxComb, ncol=length(labIds))
diffSubCt = matrix(data=NA, nrow=maxComb, ncol=length(labIds))
diffCt = matrix(data=NA, nrow=maxComb, ncol=length(labIds))
diffLabB = matrix(data=NA, nrow=maxComb, ncol=length(labIds))
diffLabSc = matrix(data=NA, nrow=maxComb, ncol=length(labIds))
diffBVsSc = matrix(data=NA, nrow=maxComb, ncol=length(labIds))
#so, these are matrices where each column represents a lab
#the matrices are just filled with values from pairs, as many as are found.

#loop through all labs
for (lab in labIds) {
  #keep track of how the next free index in each vector
  ti = 1
  si = 1
  di = 1
  tissi = 1
  subcti = 1
  cti = 1
  labbi = 1
  labsci = 1
  bvssci = 1
  #loop through all combinations
  for (samp1 in 1:(length(labs)-1)) {
    if (labs[samp1] == lab) {
      for (samp2 in (samp1+1):length(labs)) {
        #for individual and technical replicate
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
        #for tissue
        if ((labs[samp2] == lab) & 
            (cellTypes[samp1] == cellTypes[samp2]) &
            (subCellTypes[samp1] == subCellTypes[samp2]) &
            (is.na(techRepl[samp1]) | is.na(techRepl[samp2]) | (techRepl[samp1] != techRepl[samp2])) & #not technical replicates
            (is.na(individual[samp1]) | is.na(individual[samp2]) | (individual[samp1] != individual[samp2])) & #different individual
            (tissues[samp1] != tissues[samp2]) #different tissue
            ) {
          #calculate difference
          diff = sd(log2((tmmScAndBulk[,samp1] + pseudoCount)/(tmmScAndBulk[,samp2] + pseudoCount)))
          diffTissue[tissi,lab] = diff
          tissi = tissi+1
        }
        #for sub cell type
        if ((labs[samp2] == lab) & 
            (cellTypes[samp1] == cellTypes[samp2]) &
            (subCellTypes[samp1] != subCellTypes[samp2]) &
            (is.na(techRepl[samp1]) | is.na(techRepl[samp2]) | (techRepl[samp1] != techRepl[samp2])) & #not technical replicates
            (is.na(individual[samp1]) | is.na(individual[samp2]) | (individual[samp1] != individual[samp2])) & #different individual
            (tissues[samp1] == tissues[samp2])
            ) {
          #calculate difference
          diff = sd(log2((tmmScAndBulk[,samp1] + pseudoCount)/(tmmScAndBulk[,samp2] + pseudoCount)))
          diffSubCt[subcti,lab] = diff
          subcti = subcti+1
        }
        #for cell type - ignore sub cell type here
        if ((labs[samp2] == lab) & 
            (cellTypes[samp1] != cellTypes[samp2]) &
            (is.na(techRepl[samp1]) | is.na(techRepl[samp2]) | (techRepl[samp1] != techRepl[samp2])) & #not technical replicates
            (is.na(individual[samp1]) | is.na(individual[samp2]) | (individual[samp1] != individual[samp2])) & #different individual
            (tissues[samp1] == tissues[samp2])
            ) {
          #calculate difference
          diff = sd(log2((tmmScAndBulk[,samp1] + pseudoCount)/(tmmScAndBulk[,samp2] + pseudoCount)))
          diffCt[cti,lab] = diff
          cti = cti+1
        }
        #for different labs, comparison within bulk
        if ((labs[samp2] != lab) &
            ((scOrBulk[samp1] == scOrBulk[samp2]) & (scOrBulk[samp1] == 1)) & #both are bulk
            (cellTypes[samp1] == cellTypes[samp2]) &
            (subCellTypes[samp1] == subCellTypes[samp2]) &
            (is.na(techRepl[samp1]) | is.na(techRepl[samp2]) | (techRepl[samp1] != techRepl[samp2])) & #not technical replicates
            (is.na(individual[samp1]) | is.na(individual[samp2]) | (individual[samp1] != individual[samp2])) & #different individual
            (tissues[samp1] == tissues[samp2])
            ) {
          #calculate difference
          diff = sd(log2((tmmScAndBulk[,samp1] + pseudoCount)/(tmmScAndBulk[,samp2] + pseudoCount)))
          diffLabB[labbi,lab] = diff
          labbi = labbi+1
        }
        #for different labs, comparison within sc
        if ((labs[samp2] != lab) &
            ((scOrBulk[samp1] == scOrBulk[samp2]) & (scOrBulk[samp1] == 0)) & #both are sc
            (cellTypes[samp1] == cellTypes[samp2]) &
            (subCellTypes[samp1] == subCellTypes[samp2]) &
            (is.na(techRepl[samp1]) | is.na(techRepl[samp2]) | (techRepl[samp1] != techRepl[samp2])) & #not technical replicates
            (is.na(individual[samp1]) | is.na(individual[samp2]) | (individual[samp1] != individual[samp2])) & #different individual
            (tissues[samp1] == tissues[samp2])
            ) {
          #calculate difference
          diff = sd(log2((tmmScAndBulk[,samp1] + pseudoCount)/(tmmScAndBulk[,samp2] + pseudoCount)))
          diffLabSc[labsci,lab] = diff
          labsci = labsci+1
        }
        #sc vs bulk
        if ((labs[samp2] != lab) &
            (scOrBulk[samp1] != scOrBulk[samp2]) &
            (cellTypes[samp1] == cellTypes[samp2]) &
            (subCellTypes[samp1] == subCellTypes[samp2]) &
            (is.na(techRepl[samp1]) | is.na(techRepl[samp2]) | (techRepl[samp1] != techRepl[samp2])) & #not technical replicates
            (is.na(individual[samp1]) | is.na(individual[samp2]) | (individual[samp1] != individual[samp2])) & #different individual
            (tissues[samp1] == tissues[samp2])
            ) {
          #calculate difference
          diff = sd(log2((tmmScAndBulk[,samp1] + pseudoCount)/(tmmScAndBulk[,samp2] + pseudoCount)))
          diffBVsSc[bvssci,lab] = diff
          bvssci = bvssci+1
        }
      }
    }
  }
}


library(ggplot2)

##################################
#Assemble the populations to plot
##################################

techRepLab3 = tecRep[!is.na(tecRep[,3]),3]
diffLab3 = diffInd[!is.na(diffInd[,3]),3]
techRepfact = mean(techRepLab3)/mean(diffLab3)

sameIndLab4 = sameInd[!is.na(sameInd[,4]),4]
diffIndLab4 = diffInd[!is.na(diffInd[,4]),4]
sameIndfact = mean(sameIndLab4)/mean(diffIndLab4)

diffTissueLab5 = diffTissue[!is.na(diffTissue[,5]),5]
diffIndLab5 = diffInd[!is.na(diffInd[,5]),5]
diffSubCtLab5 = diffSubCt[!is.na(diffSubCt[,5]),5]
diffCtLab5 = diffCt[!is.na(diffCt[,5]),5] #There are plenty of values for other labs as well here...
diffCtAllBulk = c(diffCt[!is.na(diffCt[,1]),1],
                  diffCt[!is.na(diffCt[,2]),2],
                  diffCt[!is.na(diffCt[,4]),4],
                  diffCt[!is.na(diffCt[,5]),5])
diffCtAllSc   = c(diffCt[!is.na(diffCt[,6]),6],
                  diffCt[!is.na(diffCt[,7]),7],
                  diffCt[!is.na(diffCt[,8]),8],
                  diffCt[!is.na(diffCt[,9]),9],
                  diffCt[!is.na(diffCt[,10]),10])
diffLabAllBulk = c(diffLabB[!is.na(diffLabB[,2]),2],
                   diffLabB[!is.na(diffLabB[,3]),3],
                   diffLabB[!is.na(diffLabB[,4]),4],
                   diffLabB[!is.na(diffLabB[,5]),5])
diffLabAllSc   = c(diffLabSc[!is.na(diffLabSc[,6]),6],
                   diffLabSc[!is.na(diffLabSc[,7]),7],
                   diffLabSc[!is.na(diffLabSc[,8]),8],
                   diffLabSc[!is.na(diffLabSc[,9]),9])
diffAllBVsSc  = c(diffBVsSc[!is.na(diffBVsSc[,1]),1],
                  diffBVsSc[!is.na(diffBVsSc[,2]),2],
                  diffBVsSc[!is.na(diffBVsSc[,3]),3],
                  diffBVsSc[!is.na(diffBVsSc[,4]),4],
                  diffBVsSc[!is.na(diffBVsSc[,5]),5])
diffBulkCtrl  = c(diffInd[!is.na(diffInd[,1]),1],
                  diffInd[!is.na(diffInd[,2]),2],
                  diffInd[!is.na(diffInd[,3]),3],
                  diffInd[!is.na(diffInd[,4]),4],
                  diffInd[!is.na(diffInd[,5]),5])
diffScCtrl    = c(diffInd[!is.na(diffInd[,6]),6],
                  diffInd[!is.na(diffInd[,7]),7])



joinedData = c(techRepLab3, 
               diffLab3, 
               sameIndLab4, 
               diffIndLab4, 
               diffTissueLab5, 
               diffSubCtLab5, 
               diffIndLab5, 
               diffCtLab5,
               diffCtAllBulk,
               diffCtAllSc,
               diffLabAllBulk,
               diffLabAllSc,
               diffAllBVsSc,
               diffBulkCtrl,
               diffScCtrl)
boxes = c(rep(1,length(techRepLab3)), 
          rep(2,length(diffLab3)), 
          rep(3,length(sameIndLab4)), 
          rep(4,length(diffIndLab4)), 
          rep(5,length(diffTissueLab5)), 
          rep(6,length(diffSubCtLab5)),
          rep(7,length(diffIndLab5)), 
          rep(8,length(diffCtLab5)),
          rep(9,length(diffCtAllBulk)),
          rep(10,length(diffCtAllSc)),
          rep(11,length(diffLabAllBulk)),
          rep(12,length(diffLabAllSc)),
          rep(13,length(diffAllBVsSc)),
          rep(14,length(diffBulkCtrl)),
          rep(15,length(diffScCtrl)))

boxFac = factor(boxes, 1:15, c("Techn. repl. Bulk3", 
                               "Diff. indiv. Bulk3", 
                               "Same indiv. Bulk4", 
                               "Diff. indiv. Bulk4", 
                               "Diff. Tissue Bulk5", 
                               "Cell subtype Bulk5", 
                               "Diff. indiv. Bulk5", 
                               "Cell type Bulk5",
                               "Cell type all bulk",
                               "Cell type all sc",
                               "Lab all bulk",
                               "Lab all sc",
                               "Sc vs bulk",
                               "Bulk diff. indiv.",
                               "Sc diff. indiv."))

linFilt = boxes == 2 | boxes == 4 | boxes == 8


boxPlotWithDots <- function(data, boxes, lin){
  
  #randomize x:es:
  xes = runif(length(data)) - 0.5
  df <- data.frame(d=data, x=xes, boxes = boxes, lin=lin)
  
  g1 <- ggplot(df, aes(y=d)) + 
    geom_point(alpha=0.5, aes(x=x, color=boxes),size=2) + 
    geom_boxplot(aes(fill=NA), outlier.shape = NA) +
    scale_fill_manual(values = alpha(c("blue"), 0.0))+
    labs(y="Std(log fold change)", x="")
  g1 <- g1 + geom_segment(data=df[df$lin,], aes(x=0.7, y=0, xend=0.7, yend=3), color = "black", linetype = "dashed")

  g1 = g1+theme(panel.grid.major= element_blank(),
                panel.grid.minor= element_blank(),
                panel.background= element_blank(),
                panel.border= element_blank(),
                legend.position="none",
                axis.text.x=element_blank(), 
                axis.ticks.x=element_blank(), 
                axis.ticks.y=element_blank(),
                panel.spacing = unit(0.1, "lines"))
}


p = boxPlotWithDots(joinedData, boxFac, linFilt)
p<-p + facet_wrap( ~ boxes, nrow=1, strip.position = "bottom") +
  theme(strip.text.x = element_text(size=10, angle=90, hjust=1),
        strip.background = element_rect(colour="transparent", fill="transparent"))



annotate_figure(p,
                top = text_grob("Pairwise Estimation of Variation Factors", face = "bold", size = 14))




###########################
## Figure X Bell Curves
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

ggplot(data=df, aes(x=x, y=y, group=comparisons)) +
  geom_line(aes(colour=comparisons)) +
  labs(title="Histogram for Gene Log Fold Changes Between Samples", y="Number of genes", x="Log2 fold change")


