# This file contains code for gathering data for testing to regress out UMI copies, 
# gene length and gc content x2

dataFolder = "C:/Work/R/RNASeqCTProfileEval/"



#Download instructions for HCA single-cell evaluation:
# The raw data are in GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132044.
# The processed data are in the Single Cell Portal.
# For cortex: https://portals.broadinstitute.org/single_cell/study/SCP425/single-cell-comparison-cortex-data
# For Mixture: https://portals.broadinstitute.org/single_cell/study/SCP426/single-cell-comparison-mixture-data
# For PBMC: https://portals.broadinstitute.org/single_cell/study/SCP424/single-cell-comparison-pbmc-data
# You will need cortex only. Unzip everything. Then,
# copy some files:
# genes.txt -> counts/genes.tsv, genes.txt -> UMIs/genes.tsv  - According to instructions, we should have picked genes.count.txt, but it is identical but weird it seems
# counts.reads.txt -> counts/matrix.mtx, counts.umis.txt -> UMIs/matrix.mtx
# cell.names.txt -> counts/barcodes.tsv, cell.names.txt -> UMIs/barcodes.tsv



##################################
## The goal is to create a file per sample.
## The rows are genes, not transcripts, so gene length etc will be an average
## over all transcripts for the gene
## Do this for cortex 1 and 2, none of the other have bulk data as reference.
##################################

#################
#load cortex data
#################

HCASCE_counts <- Read10X(data.dir = "C:/Work/MatlabCode/components/SingleCellToolbox/ImportableData/HCA_single-cell_comparison/cortex/counts")
HCASCE_UMIs <- Read10X(data.dir = "C:/Work/MatlabCode/components/SingleCellToolbox/ImportableData/HCA_single-cell_comparison/cortex/UMIs")

HCASCE_bulk1 = read.table("C:/Work/MatlabCode/components/SingleCellToolbox/ImportableData/HCA_single-cell_comparison/cortex/bulk/cortex1.genes.results",header=T, sep="\t", row.names = 1)
HCASCE_bulk2 = read.table("C:/Work/MatlabCode/components/SingleCellToolbox/ImportableData/HCA_single-cell_comparison/cortex/bulk/cortex2.genes.results",header=T, sep="\t", row.names = 1)

#convert gene ids to ensembl, this simplifies further down
library(biomaRt)
ensembl_us_west = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="uswest.ensembl.org")
listDatasets(ensembl_us_west)
ensembl = useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")

geneConvTableM <- getBM(attributes=c('ensembl_gene_id', 'mgi_symbol'), mart = ensembl)


umiGenes = row.names(HCASCE_UMIs)
ind = match(umiGenes,geneConvTableM$mgi_symbol)
newUMIGenes = geneConvTableM$ensembl_gene_id[ind]


cortexBulkTPM = cbind(HCASCE_bulk1$TPM,HCASCE_bulk2$TPM)
colnames(cortexBulkTPM) = cbind("Bulk1", "Bulk2")
row.names(cortexBulkTPM) = row.names(HCASCE_bulk1)


#now, handle 10x: UMIs and counts compared to bulk TPM
cortex1UMIs = as.matrix(HCASCE_UMIs[,grep("Cortex1_10xChromium", colnames(HCASCE_UMIs))])
cortex2UMIs = as.matrix(HCASCE_UMIs[,grep("Cortex2_10xChromium", colnames(HCASCE_UMIs))])
cortex1Counts = as.matrix(HCASCE_counts[,grep("Cortex1_10xChromium", colnames(HCASCE_counts))])
cortex2Counts = as.matrix(HCASCE_counts[,grep("Cortex2_10xChromium", colnames(HCASCE_counts))])
cortex12UMIAndCounts = cbind(rowSums(cortex1UMIs), rowSums(cortex2UMIs), rowSums(cortex1Counts), rowSums(cortex2Counts))
colnames(cortex12UMIAndCounts) = c("UMI1","UMI2", "count1", "count2")
row.names(cortex12UMIAndCounts) = newUMIGenes;
#remove genes for which conversion failed: 
cortex12UMIAndCounts = cortex12UMIAndCounts[!is.na(row.names(cortex12UMIAndCounts)),]

rm(cortex1UMIs, cortex2UMIs, cortex1Counts, cortex2Counts)

#merge UMI and bulk - don't use merge, it is super slow
library(dplyr)
library(tibble)
cortex12Merged = inner_join(rownames_to_column(as.data.frame(cortexBulkTPM)), rownames_to_column(as.data.frame(cortex12UMIAndCounts)))
row.names(cortex12Merged) = cortex12Merged$rowname
cortex12Merged = cortex12Merged[,-1]
#cortex12MergedTPM = MakeTPM(cortex12Merged)
#cortex12MergedTMM = TMMNorm(cortex12MergedTPM)


###############
## Now get UMI copies, gene length and GC content
###############



###############
# Gene length
##############
#BiocManager::install(c("GenomicFeatures"))
library("GenomicFeatures")

#txdb = makeTxDbFromBiomart(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", transcript_ids=NULL, circ_seqs=DEFAULT_CIRC_SEQS, filter=NULL, id_prefix="ensembl_", host="www.ensembl.org", port=80, taxonomyId=NA, miRBaseBuild=NA)
#GRCm38.p6
txdbM = makeTxDbFromBiomart(biomart="ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl", transcript_ids=NULL, circ_seqs=DEFAULT_CIRC_SEQS, filter=NULL, id_prefix="ensembl_", host="www.ensembl.org", port=80, taxonomyId=NA, miRBaseBuild=NA)

tlM = transcriptLengths(txdbM, with.cds_len=FALSE)
row.names(tlM) = tlM$tx_name


##############
#GC Content
##############

#load genome for mus musculus
#BiocManager::install("BSgenome")
#BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")
library(BSgenome.Mmusculus.UCSC.mm10)
genomeM <- BSgenome.Mmusculus.UCSC.mm10
BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdbM <- TxDb.Mmusculus.UCSC.mm10.knownGene
transcriptsM <- exonsBy(txdbM, by="tx", use.names=TRUE)


seqsM = extractTranscriptSeqs(genomeM, transcriptsM)

#calculate gc content in the strings

gcContent <- function(x) {
  return (letterFrequency(x, c("GC"), OR="|", as.prob=TRUE))
} 

gcContentTail <- function(sequences, len) {
  #so, first reverse, then loop through all and
  #pick the first up to length. If shorter, use the
  #shorter value. Had to make a loop, couldn't make vapply work!
  tmp = reverse(sequences)
  #couldn't make it work with vapply, using a loop
  #this takes 45 min to run
  for (i in 1:length(sequences)) {
    d = min(width(tmp)[[i]],len)
    tmp[i] = subseq(tmp[i],1,d)
    if (i %% 1000 == 0) {
      print(100*i/length(sequences))
    }
  }
  
  return (gcContent(tmp))
}



#test
#gcContent
b = BString("GGCCGA")
gcContent(b)#should be 5/6 = 0.8333333, ok!
#gcContentTail
#the last 15 letters of the first gene is: "ACCTTTGCATATAAA", so this should be 4/15 = 0.2666667
test = gcContentTail(seqsM[1],15) #ok


gcFullLength = gcContent(seqsM)
gcTail = gcContentTail(seqsM, 150)

#these have transcript ids only, not gene ids. So, merge with the length to get gene id
gcs = data.frame(gcFullLength, gcTail)
colnames(gcs) = c("gcFullLength", "gcTail")
genes = seqsM@ranges@NAMES
#need to remove the version from the transcript name
genes = substr(genes, 1, 18)

row.names(gcs) = genes

gcsMMerged = inner_join(rownames_to_column(gcs), rownames_to_column(tlM))
row.names(gcsMMerged) = gcsMMerged$rowname
gcsMMerged = gcsMMerged[,-1]

#take mean of all transcripts for each gene
gcsM = aggregate(gcsMMerged[,c(1,2,7)], list(gcsMMerged[,5]), mean)
row.names(gcsM) = gcsM$Group.1
gcsM = gcsM[,-1]

#now merge with the data
cortex12TotMerged = inner_join(rownames_to_column(as.data.frame(cortex12Merged)), rownames_to_column(as.data.frame(gcsM)))
row.names(cortex12TotMerged) = cortex12TotMerged$rowname
cortex12TotMerged = cortex12TotMerged[,-1]

#assumes the following structure: Bulk1 Bulk2 UMI1 UMI2 count1 count2 gcFullLength    gcTail   tx_len
extractSample <- function(mergedData, index) {
  addN = index - 1
  dat = mergedData[,c(1+addN, 3+addN, 7, 8, 9)]
  #calculate removed counts' fraction
  dat = cbind(dat, (mergedData[,5+addN]- mergedData[,3+addN]) / mergedData[,5+addN])
  #add TPM, and TMM-normalized, log transformed data
  d2 = dat[,c(1,2)]
  d2 = MakeTPM(d2);
  tmmNorm = TMMNorm(d2)
  d3 = log2(tmmNorm + 0.05)
  d3 = cbind(d3, log2((tmmNorm[,2] + 0.05)/(tmmNorm[,1] + 0.05)))
  dat = cbind(dat, d2,d3)
  colnames(dat) = c("bulk", "UMI", "gcFullLength", "gcTail", "geneLength", "remUMIFrac", "bulkTPM", "UMITPM", "logBulkTMM", "logUMITMM", "LogUMIDivBulk")
  return (dat)
}

cortex1Data = extractSample(cortex12TotMerged, 1)
cortex2Data = extractSample(cortex12TotMerged, 2)

#save to disk
write.table(cortex1Data, file=paste0(dataFolder, "/ScVsBulkCortex1.txt"), row.names=T, sep="\t")
write.table(cortex2Data, file=paste0(dataFolder, "/ScVsBulkCortex2.txt"), row.names=T, sep="\t")






