install.packages("devtools")
install.packages("knitr")
devtools::install_github("GfellerLab/EPIC", build_vignettes=TRUE)
library("EPIC")

#loop through all LUSC TCGA files
#folder = "C:/Work/MatlabCode/projects/AgileCancerTreatment/data/TCGAFiles/TCGA-LUSC"

#out.file<-""

#file.names <- dir(folder, pattern =".txt")
#for(i in 1:length(file.names)){
#for(i in 1){
#  data <- read.table(paste0(folder, "/", file.names[i]),header=FALSE, sep="\t", stringsAsFactors=FALSE)
  #out.file <- rbind(out.file, file)
#}


#read the file where all the samples are merged into one and converted to HGNC (in the cpp code, we
#may need to fix this later)

data = read.csv("C:/Work/MatlabCode/projects/AgileCancerTreatment/data/TCGAFiles/TCGA_SKCM.txt", header=TRUE, row.names=1, sep = "\t", quote = "\"'", check.names=FALSE);

print('data reading done')

print (data[4,1])#expected 0.2852301

#use row.names(data) to get the genes and colnames(data) to get the sample ids
mat <- data.matrix(data)
deconvResults <- EPIC(bulk = mat)
#deconvResults <- EPIC(bulk = mat);
#write.matrix(deconvResults$mRNAProportions, file = "DeconvResults_SKCM.txt", sep = "\t");
write.table(deconvResults$mRNAProportions, file = "../../data/TCGAFiles/DeconvResults_SKCM.txt", sep = "\t", col.names = NA, quote=FALSE)

#use the line below to write the signature genes to file
#write(TRef$sigGenes, file = "C:/Work/MatlabCode/projects/AgileCancerTreatment/data/EPIC/SigGenes.txt", ncolumns = 1, append = FALSE, sep = "\t");

head(deconvResults$mRNAProportions)
colMeans(deconvResults$mRNAProportions)
head(mat)
rowSums(mat)

covariates = deconvResults$mRNAProportions;

#Loop through all genes. For each gene, make a linear interpolation with the 
#fraction of different covariates 
#for (i in 1:(dim(mat)[1])) {
for (i in 23001) {
  expr = mat[i,]
  expr = mat[row.names(mat) == "CD8A",]
  res = lm(expr ~ covariates)
}

#try non-negative least squares:
#install.packages("nnls")
#library("nnls")
res2 = nnls(covariates, expr)
res2



