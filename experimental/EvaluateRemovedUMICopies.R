
library("ggplot2")
library("Seurat")


cts <- Read10X(data.dir = "C:/Work/MatlabCode/components/SingleCellToolbox/ImportableData/HCA_single-cell_comparison/cortex/counts")
umis <- Read10X(data.dir = "C:/Work/MatlabCode/components/SingleCellToolbox/ImportableData/HCA_single-cell_comparison/cortex/UMIs")

cts10xc1 = as.matrix(cts[,grep("Cortex1_10xChromium", colnames(cts))])
umis10xc1 = as.matrix(umis[,grep("Cortex1_10xChromium", colnames(umis))])
diff = cts10xc1 - umis10xc1

sel = cts10xc1 > 0

c1counts = cts10xc1[sel]
c1umis = umis10xc1[sel]
c1diff = diff[sel]

countsPerUmi = c1counts / c1umis
hist(countsPerUmi, (0:45)-0.5)

max(countsPerUmi)

sum(c1umis[c1umis>1])
#5370376
sum(c1umis[c1umis==1])
#2380634


sel2 = umis10xc1 == 1
c1counts2 = cts10xc1[sel2]
c1umis2 = umis10xc1[sel2]
c1diff2 = diff[sel2]
countsPerUmi2 = c1counts2 / c1umis2
hist(countsPerUmi2, (0:45)-0.5)


#so, it seems like 30.7 % of the reads are ones, that cannot have a duplicated UMI due to 
#read errors