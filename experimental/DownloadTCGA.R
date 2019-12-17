library("data.table")
library("readr")
library("TCGAbiolinks")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("TCGAbiolinks")


query <- GDCquery(project = "TCGA-LUAD", sample.type="Primary solid Tumor",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - FPKM")

GDCdownload(query, directory = "C:/Work/R/RNASeqCTProfileEval/TCGA/LUAD")

res = query$results[[1]]

fns = res$file_name

endChars = regexpr("\\.", fns) - 1
fileIds = mapply(substring, fns, 1, endChars)


id2BarcodeMapping = cbind(fileIds, res$cases)


#after this is done:
#1. Use search in the explorer on *.gz files, mark them all
#2. Use 7zip to extract, select extract files, change the path to
#   [path to repo]/TCGA/LUAD/ and untick the box below the path
#3. Run the matlab script AssembleTCGALUAD.m

#now download purity

#first get the first line from the data table
assembledData = read.table("C:/Work/R/RNASeqCTProfileEval/TCGA/LUAD.txt")


ids = colnames(assembledData)
ids = gsub("\\.", "-", ids) # '-' changed to '.'
ids = gsub("X", "", ids) # X added for names starting with a number


ind = match(ids,id2BarcodeMapping[,1])
barcodes = id2BarcodeMapping[ind, 2] #the barcodes here are in the same order as the samples
barcodes = mapply(substring, barcodes, 1, 16)#only take the first 16 chars of the barcode



purind = match(barcodes,Tumor.purity[,1])
purity = Tumor.purity[purind,]
purity

#load the CIBERSORT results and show correlation with purity

deconvResults = read.table("C:/Work/R/RNASeqCTProfileEval/CIBERSORTx_Job9_Adjusted.txt", header=T, row.names = 1)


ind = match(deconvResults$Mixture,id2BarcodeMapping[,1])
barcodes = id2BarcodeMapping[ind, 2] #the barcodes here are in the same order as the samples
barcodes = mapply(substring, barcodes, 1, 16)#only take the first 16 chars of the barcode


purind = match(barcodes,Tumor.purity[,1])
purity = Tumor.purity[purind,]

mal = as.numeric(deconvResults$Malignant)
#for some reason the purity is a factor, with decimal numbers using ",", probably some locale messing it up
#so, they need to be converted
pur = as.character(purity$CPE)
pur = gsub(",", ".", pur) # replace , with .
pur2 = as.numeric(pur)
plot(pur2, mal)

#test to run epic on the samples
#...



