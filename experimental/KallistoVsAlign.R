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
