# This file contains code for gathering data for testing to regress out UMI copies, 
# gene length and gc content x2

dataFolder = "C:/Work/R/RNASeqCTProfileEval/OrigBulkData/"

################
# BLUEPRINT sample
################

bp = read.table(paste0(dataFolder, "C001FRB1.gene_quantification.rsem_grape2_crg.GRCh38.20150622.results"), header = T,sep="\t")
bpgenes = substr(bp[,1], 1, 15) #remove version on the genes
bpexpr = bp[,6]

#convert to gene names
library(biomaRt)
ensembl_us_west = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="uswest.ensembl.org")
#listDatasets(ensembl_us_west)
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
geneConvTableH <- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol'), mart = ensembl)

ind = match(bpgenes, geneConvTableH$ensembl_gene_id)
newGenes = geneConvTableH$hgnc_symbol[ind]
sel = (newGenes != "" & !is.na(newGenes) ) & !duplicated(newGenes) # also remove duplicated genes, difficult to know how to interpret those
bpexpr = bpexpr[sel]
bpgenes = newGenes[sel]

bpdata = data.frame(bpexpr)
row.names(bpdata) = bpgenes

################
# ENCODE sample
################

enc = read.table(paste0(dataFolder, "ENCFF088DIY.tsv"), header = T,sep="\t")
encgenes = substr(enc[,1], 1, 15) #remove version on the genes
encexpr = enc[,6]

ind = match(encgenes, geneConvTableH$ensembl_gene_id)
newGenes = geneConvTableH$hgnc_symbol[ind]
sel = (newGenes != "" & !is.na(newGenes) ) & !duplicated(newGenes)
encexpr = encexpr[sel]
encgenes = newGenes[sel]

encdata = data.frame(encexpr)
row.names(encdata) = encgenes


################
# FANTOM5 sample
################

fantom5 = read.table(paste0(dataFolder, "hg19.gene_phase1and2combined_tpm.osc.txt"), header = T,sep="\t", comment.char = "#", as.is=c(1))

fantom5genes = fantom5[,1]
fantom5expr = fantom5$tpm.CD8.2b.20T.20Cells.20.28pluriselect.29.2c.20donor090309.2c.20donation1.CNhs12176.12186.129A8

fantom5data = data.frame(fantom5expr)
row.names(fantom5data) = fantom5genes


#######################
# Now merge the data
#######################

interMed = inner_join(rownames_to_column(as.data.frame(bpdata)), rownames_to_column(as.data.frame(encdata)))
row.names(interMed) = interMed$rowname
interMed = interMed[,-1]

joinedData = inner_join(rownames_to_column(as.data.frame(interMed)), rownames_to_column(as.data.frame(fantom5data)))
row.names(joinedData) = joinedData$rowname
joinedData = joinedData[,-1]



joinedData = MakeTPM(joinedData)

write.table(joinedData, file=paste0(dataFolder, "../AuthorsProcBulk.txt"), row.names=T, sep="\t")




