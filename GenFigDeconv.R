
library(ggplot2)
library(readr)
fig_path = "Z:/projects/Cell type profiles/figures/"


#read data
#########################

labRes = c("l1", "l4", "l5", "l6", "l7")
scRes = c("pbmc68k", "mel", "lct", "lch", "hca")
methods = c("", "Q", "BatchB", "BatchS")

numRuns = 5*3 + 5*4
sampPerRun = 40

resultsPath = "C:/Work/R/RNASeqCTProfileEval/data/deconv/results/"

results = matrix(0, nrow=sampPerRun, ncol=numRuns)


for (i in 1:length(labRes)) {
  for (m in 1:3) {
    #load the results file
    f = read_tsv(file=paste0(resultsPath, labRes[i], methods[m], ".txt"))
    results[,(i-1)*3+m] = f$BCells
  }
}

for (i in 1:length(scRes)) {
  for (m in 1:4) {
    #load the results file
    f = read_tsv(file=paste0(resultsPath, scRes[i], methods[m], ".txt"))
    results[,15+(i-1)*4+m] = f$BCell
  }
}

resultsInt = matrix(0, nrow=9, ncol=3)
resultsInt5 = matrix(0, nrow=8, ncol=3)


for (m in 1:3) {
  #load the results file
  f = read_tsv(file=paste0(resultsPath, "l4Int", methods[m], ".txt"))
  resultsInt[,m] = f$BCells
}

for (m in 1:3) {
  #load the results file
  f = read_tsv(file=paste0(resultsPath, "l5Int", methods[m], ".txt"))
  resultsInt5[,m] = f$BCells
}

resultsLM22 = matrix(0, nrow=sampPerRun, ncol=3)
for (m in 1:3) {
  #load the results file
  f = read_tsv(file=paste0(resultsPath, "lm22", methods[m], ".txt"))
  bcells = f$`B cells naive` + f$`B cells memory` + f$`Plasma cells`
  
  resultsLM22[,m] = bcells
}


# Now plot the results
#######################

#Use a grouped boxplot

# create a data frame
profIds = c(rep(4:8, each=3), rep(9:13, each=4))
meth = c(rep(1:3, 5), rep(1:4, 5))

#generate a melted data frame
df = data.frame(bfrac = rep(0, sampPerRun * numRuns + 9*3 + 8*3 + sampPerRun * 3), 
                profIds = c(rep(profIds, each=sampPerRun), rep(2,9*3), rep(3,8*3), rep(1,sampPerRun*3)), 
                meth = c(rep(meth, each=sampPerRun), rep(1:3, each=9), rep(1:3, each=8), rep(1:3, each=sampPerRun)))
#fill in the values
for (col in 1:numRuns) {
  for (row in 1:sampPerRun) {
    df$bfrac[(col-1)*sampPerRun + row] = results[row, col]
  }  
}

#add the internal values
for (col in 1:3) {
  for (row in 1:9) {
    df$bfrac[numRuns*sampPerRun + (col-1)*9 + row] = resultsInt[row, col]
  }  
}

#add the internal values lab 5
for (col in 1:3) {
  for (row in 1:8) {
    df$bfrac[numRuns*sampPerRun + 3*9 + (col-1)*8 + row] = resultsInt5[row, col]
  }  
}

#and LM22
for (col in 1:3) {
  for (row in 1:sampPerRun) {
    df$bfrac[numRuns*sampPerRun + 3*9 + 3*8 + (col-1)*sampPerRun + row] = resultsLM22[row, col]
  }  
}


df$profIds = factor(df$profIds, 1:13, 
                    c("1: LM22", "2: Bulk 4 Internal","3: Bulk 5 Internal","4: Bulk 1", "5: Bulk 4", "6: Bulk 5", "7: Pooled SC HCA CB", "8: Pooled SC LC", 
                      "9: SC PBMC68k", "10: SC Melanoma", "11: SC LC Tumor", "12: SC LC Healthy Tis.", "13: SC HCA CB"))
df$meth = factor(df$meth, 1:4, 
                 c("TPM/CPM", "Quantile", "Batch Corr B", "Batch Corr S"))



#now filter out the results that belong to the same lab as the profiles are generated from
#Bulk 1 (sample 1-5)
toRemBulk1 = c(1:5, 1:5 + sampPerRun, 1:5 + 2 * sampPerRun)
toRemBulk4 = c(12:17 + 3*sampPerRun, 12:17 + 4*sampPerRun, 12:17 + 5*sampPerRun)
toRemBulk5 = c(18:40 + sampPerRun*6, 18:40 + sampPerRun*7, 18:40 + sampPerRun*8)
toRem = c(toRemBulk1, toRemBulk4, toRemBulk5)
dfFilt = df[-toRem,]


#recalculate data to relative error
dfFilt$bfrac = abs(dfFilt$bfrac - 0.5) / 0.5

#check how many mixtures we have for each run (multiplied with the number of methods)
dfFilt %>% group_by(profIds) %>% tally()

# grouped boxplot
pDeconv = ggplot(dfFilt, aes(x=profIds, y=bfrac, fill=meth)) + 
#  geom_boxplot(outlier.shape = 21) +
  geom_boxplot(outlier.shape = NA, coef = 0) +
  ylab("Relative Error") + xlab('Cell-type profile set') + 
  #ggtitle("Deconvolution Performance") +
  theme(axis.text.x = element_text( #size  = 11,
                                   angle = 20,
                                   hjust = 1,
                                   vjust = 1), 
        legend.title = element_blank(), 
        legend.position="bottom",
        panel.background = element_rect("white", "white", 0, 0, "white"))


#skip title
fig7 = pDeconv
#fig7 = annotate_figure(pDeconv,
#                       top = text_grob("Deconvolution Performance", face = "bold", size = 14))

fig7


ggsave(
  paste0(fig_path, "Fig7Plot.png"),
  plot = fig7,
  width = 6, height = 4, dpi = 300)


#supporting with outliers and whiskers
pDeconvS = ggplot(dfFilt, aes(x=profIds, y=bfrac, fill=meth)) + 
  geom_boxplot(outlier.shape = 21) +
  ylab("Relative Error") + xlab('Cell-type profile set') + 
  #ggtitle("Deconvolution Performance") +
  theme(axis.text.x = element_text( #size  = 11,
    angle = 20,
    hjust = 1,
    vjust = 1), 
    legend.title = element_blank(), 
    legend.position="bottom",
    panel.background = element_rect("white", "white", 0, 0, "white"))

#skip title
figS7 = pDeconvS

figS7


ggsave(
  paste0(fig_path, "FigS7.png"),
  plot = figS7,
  width = 6, height = 4, dpi = 300)



