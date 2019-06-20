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


