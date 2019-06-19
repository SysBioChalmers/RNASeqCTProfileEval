library(plyr)

dataFolder = "E:/BulkProfiles/BulkProfiles - Copy"


#So, the files were processed in two ways:
#1. ENCODE and FAMTOM 5 samples were first copied to file names according to their 
#   place in the design matrix, and processed accordingly. Kallisto was called manually.
#2. To save disk space, this was not done for the BLUEPRINT files; they were instead processed
#   in their original structure. Kallisto was run using the python script RunKallisto.py

#process 

#for (i in 1:25) {
for (i in 25) {
  ## Read in abundance file from kallisto
  abundance = read.table(file=paste(dataFolder, "/output_",i, "/abundance.tsv", sep=""), header=T)
  head(abundance)

  # Import ENSG gene list
  gene_names = read.table(file="E:/BulkProfiles - Copy/mart_export.txt", header=T, sep="\t")
  colnames(gene_names) = c("ENSG", "HGNC", "ENST")
  colnames(abundance) = c("ENST", "length", "eff_length", "est_counts", "tpm")
  head(abundance)
  head(gene_names)

  # Align
  aligned = merge(abundance, gene_names, by = "ENST")
  head(aligned)

  # Get rid of everything but est_counts
  aligned = aligned[ ,c("ENST", "ENSG", "HGNC", "est_counts")]
  head(aligned)

  # Aggregate
  agg = aggregate(est_counts ~ HGNC, data = aligned, sum)
  head(agg)

  # Write agg to file
  write.table(agg, file=paste("E:/BulkProfiles - Copy/output_", i, "/abundance_hgnc.txt", sep=""),
              row.names=F, sep="\t")
}

# Make a big data frame with all (ENCODE and FANTOM5) samples. HGNC gene as rows and subject/sample as column

# Read in the first sample
count_df = read.table(file="E:/BulkProfiles - Copy/output_1/abundance_hgnc.txt", header=T, sep="\t")
head(count_df)

# Read in the rest
for ( i in 2:24) {
  count_df = cbind(count_df, read.table(file=paste("E:/BulkProfiles - Copy/output_",i, "/abundance_hgnc.txt", sep=""), header=T, sep="\t")[ ,"est_counts"] )
}
head(count_df)
dim(count_df)

colnames(count_df) = c("hgnc", seq(1,24,1) )
head(count_df)

# Set hgnc as rownames
rownms = count_df[ ,1]
row.names(count_df) = rownms
head(count_df)
count_df = count_df[ ,-1]
head(count_df)
dim(count_df)




############################################# READ IN Blueprint data
Blueprint_C001NBB1 = read.table(file="E:/BulkProfiles - Copy/Blueprint/BCells/EGAD00001001145/C001NBB1.gene_quantification.rsem_grape2_crg.GRCh38.20150622.results",
                                sep="\t", header=T)
head(Blueprint_C001NBB1)
dim(Blueprint_C001NBB1)
Blueprint_C001NBB1 = Blueprint_C001NBB1[ , c(1,5)] ## Keep the ENSG gene columns and the expected count
head(Blueprint_C001NBB1)
colnames(Blueprint_C001NBB1) = c("ENSG", "expected_count")
head(Blueprint_C001NBB1)

######## Convert gene id to HGNC
# Import ENSG gene list
#gene_names = read.table(file="E:/BulkProfiles - Copy/mart_export.txt", header=T, sep="\t")
#colnames(gene_names) = c("ENSG", "HGNC", "ENST")
#head(gene_names)

# Remove the ".XX" ending (IS THIS OK TO DO?)
is.data.frame(Blueprint_C001NBB1)
rownames = Blueprint_C001NBB1[ ,"ENSG"]
rownames
length(rownames)
rownames = as.vector(rownames)
is.vector(rownames)
rownames = sub('\\..*', '', rownames)
rownames
Blueprint_C001NBB1 = cbind(rownames, Blueprint_C001NBB1)
head(Blueprint_C001NBB1)
Blueprint_C001NBB1 = Blueprint_C001NBB1[ , c("rownames","expected_count")]
colnames(Blueprint_C001NBB1) = c("ENSG", "expected_count")
head(Blueprint_C001NBB1)
dim(Blueprint_C001NBB1)

# Align
aligned_blueprint = merge(Blueprint_C001NBB1, gene_names, by = "ENSG")
head(aligned_blueprint)

# Aggregate
agg_blueprint = aggregate(expected_count ~ HGNC, data = aligned_blueprint, sum)
head(agg_blueprint)
dim(agg_blueprint)






## Construct meta-data frame
Immune_cell_type = c( rep("Bcell",7), rep("Tcell",5), rep("Bcell",6), rep("Tcell",6)  )
Data_base = c( rep("ENCODE",12), rep("FANTOM5",12) )
read_type = c( rep("Paired_end",12), rep("Single_end",12) )

meta_df = cbind(Immune_cell_type, Data_base, read_type)
head(meta_df)
row.names(meta_df) = colnames(count_df)
head(meta_df)

######### RUN DESeq2

#library("DESeq2")
# Round count_df to integers
#head(count_df)
#count_df = round(count_df, 0)
#head(count_df)

#dds <- DESeqDataSetFromMatrix(countData = count_df, colData = meta_df, design = ~ Immune_cell_type)
#dds

######### RUN DESeq (because fitNbinomGLMs doesn't exist in DESeq2!)
#library(DESeq)
cds = newCountDataSet( count_df, meta_df )
cds

## Normalisation
cds = estimateSizeFactors( cds )
sizeFactors( cds )

## Variance estimation
cds = estimateDispersions( cds )

## Fit A Generalized Linear Model (GLM) For Each Gene
fit = fitNbinomGLMs( cds, count ~ Immune_cell_type*Data_base ) ## Full model
fit_reduced = fitNbinomGLMs( cds, count ~ Data_base ) ## Reduced model
p_values = nbinomGLMTest(fit, fit_reduced)

sort(p_values, decreasing=FALSE)
par(mar=c(4,4,4,4))
hist(p_values)


