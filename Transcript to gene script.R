library(plyr)



#So, the files were processed in two ways:
#1. ENCODE and FANTOM 5 samples were first copied to file names according to their 
#   place in the design matrix, i.e. 1, 2, etc, and processed accordingly. Kallisto was 
#   called manually as follows:
#   kallisto quant -i transcripts.gtf.gz -o output_4 -b 1 4_1.fastq.gz 4_2.fastq.gz
#2. To save disk space, this was not done for the BLUEPRINT files; they were instead processed
#   in their original structure. Kallisto was run using the python script RunKallisto.py

#process ENCODE and FANTOM5 

dataFolderEF = "E:/BulkProfiles/BulkProfiles - Copy"


# Import ENSG gene list
gene_names = read.table(file=paste0(dataFolderEF, "/mart_export.txt"), header=T, sep="\t")
colnames(gene_names) = c("ENSG", "HGNC", "ENST")

for (i in 1:25) {
#for (i in 1) {
  ## Read in abundance file from kallisto
  abundance = read.table(file=paste(dataFolderEF, "/output_",i, "/abundance.tsv", sep=""), header=T)
  head(abundance)

  colnames(abundance) = c("ENST", "length", "eff_length", "est_counts", "tpm")
  head(abundance)
  head(gene_names)

  # Align
  aligned = merge(abundance, gene_names, by = "ENST")
  head(aligned)

  # Get rid of everything but est_counts and tpm
  aligned = aligned[ ,c("ENST", "ENSG", "HGNC", "est_counts", "tpm")]
  head(aligned)

  # Aggregate
  agg = aggregate(. ~ HGNC, data = aligned, sum)
  agg = agg[,c('HGNC','est_counts','tpm')]
  head(agg)

  # Write agg to file
  write.table(agg, file=paste(dataFolderEF, "/output_", i, "/abundance_hgnc.txt", sep=""),
              row.names=F, sep="\t")
}

#now process the BLUEPRINT data

BPFileInfo = list(
  list("EGAD00001001137", list(
    list("EGAR00001137340_130919_SN546_0217_B_C2GFPACXX_CTTGTA_1_read1.fastq.gz", "EGAR00001137340_130919_SN546_0217_B_C2GFPACXX_CTTGTA_1_read2.fastq.gz", "EGAR00001137340_130919_SN546_0217_B_C2GFPACXX_CTTGTA_1", 26),
    list("EGAR00001219222_140618_SN935_0195_A_C42K0ACXX_ACAGTG_2_read1.fastq.gz", "EGAR00001219222_140618_SN935_0195_A_C42K0ACXX_ACAGTG_2_read2.fastq.gz", "EGAR00001219222_140618_SN935_0195_A_C42K0ACXX_ACAGTG_2", 27)
  )
  ),
  list("EGAD00001001145", list(
    list("EGAR00001074533_s_120913_6.ATTCCT.read1.fastq.gz", "EGAR00001074533_s_120913_6.ATTCCT.read2.fastq.gz", "EGAR00001074533_s_120913_6.ATTCCT", 28),
    list("EGAR00001094082_130214_SN935_0147_B_D1T9CACXX_AGTCAA_7_read1.fastq.gz", "EGAR00001094082_130214_SN935_0147_B_D1T9CACXX_AGTCAA_7_read2.fastq.gz", "EGAR00001094082_130214_SN935_0147_B_D1T9CACXX_AGTCAA_7", 29)
  )
  ),
  list("EGAD00001001173", list(
    list("EGAR00001219220_140618_SN935_0195_A_C42K0ACXX_ACTGAT_3_read1.fastq.gz", "EGAR00001219220_140618_SN935_0195_A_C42K0ACXX_ACTGAT_3_read2.fastq.gz", "EGAR00001219220_140618_SN935_0195_A_C42K0ACXX_ACTGAT_3", 30),
    list("EGAR00001246138_140814_SN546_0249_A_C48WMACXX_ACTGAT_3_read1.fastq.gz", "EGAR00001246138_140814_SN546_0249_A_C48WMACXX_ACTGAT_3_read2.fastq.gz", "EGAR00001246138_140814_SN546_0249_A_C48WMACXX_ACTGAT_3", 31),
    list("EGAR00001246140_140814_SN546_0249_A_C48WMACXX_ACAGTG_2_read1.fastq.gz", "EGAR00001246140_140814_SN546_0249_A_C48WMACXX_ACAGTG_2_read2.fastq.gz", "EGAR00001246140_140814_SN546_0249_A_C48WMACXX_ACAGTG_2", 32),
    list("EGAR00001074535_s_120913_5.GCCAAT.read1.fastq.gz", "EGAR00001074535_s_120913_5.GCCAAT.read2.fastq.gz", "EGAR00001074535_s_120913_5.GCCAAT", 33),
    list("EGAR00001219216_140618_SN935_0195_A_C42K0ACXX_GATCAG_5_read1.fastq.gz", "EGAR00001219216_140618_SN935_0195_A_C42K0ACXX_GATCAG_5_read2.fastq.gz", "EGAR00001219216_140618_SN935_0195_A_C42K0ACXX_GATCAG_5", 34),
    list("EGAR00001246139_140814_SN546_0249_A_C48WMACXX_GTGAAA_2_read1.fastq.gz", "EGAR00001246139_140814_SN546_0249_A_C48WMACXX_GTGAAA_2_read2.fastq.gz", "EGAR00001246139_140814_SN546_0249_A_C48WMACXX_GTGAAA_2", 35),
    list("EGAR00001219218_140618_SN935_0195_A_C42K0ACXX_ATCACG_4_read1.fastq.gz", "EGAR00001219218_140618_SN935_0195_A_C42K0ACXX_ATCACG_4_read2.fastq.gz", "EGAR00001219218_140618_SN935_0195_A_C42K0ACXX_ATCACG_4", 36),
    list("EGAR00001246137_140814_SN546_0249_A_C48WMACXX_GTCCGC_3_read1.fastq.gz", "EGAR00001246137_140814_SN546_0249_A_C48WMACXX_GTCCGC_3_read2.fastq.gz", "EGAR00001246137_140814_SN546_0249_A_C48WMACXX_GTCCGC_3", 37),
    list("EGAR00001137343_130919_SN546_0217_B_C2GFPACXX_GCCAAT_1_read1.fastq.gz", "EGAR00001137343_130919_SN546_0217_B_C2GFPACXX_GCCAAT_1_read2.fastq.gz", "EGAR00001137343_130919_SN546_0217_B_C2GFPACXX_GCCAAT_1", 38),
    list("EGAR00001219224_140618_SN935_0195_A_C42K0ACXX_CTTGTA_1_read1.fastq.gz", "EGAR00001219224_140618_SN935_0195_A_C42K0ACXX_CTTGTA_1_read2.fastq.gz", "EGAR00001219224_140618_SN935_0195_A_C42K0ACXX_CTTGTA_1", 39)
  )
  ),
  list("EGAD00001001478", list(
    list("EGAR00001270597_141030_SN935_0211_B_C5R24ACXX_GTCCGC_3_read1.fastq.gz", "EGAR00001270597_141030_SN935_0211_B_C5R24ACXX_GTCCGC_3_read2.fastq.gz", "EGAR00001270597_141030_SN935_0211_B_C5R24ACXX_GTCCGC_3", 40)
  )
  ),
  list("EGAD00001001483", list(
    list("EGAR00001270597_141030_SN935_0211_B_C5R24ACXX_GTCCGC_3_read1.fastq.gz", "EGAR00001270597_141030_SN935_0211_B_C5R24ACXX_GTCCGC_3_read2.fastq.gz", "EGAR00001270597_141030_SN935_0211_B_C5R24ACXX_GTCCGC_3", 41)
  )
  ),
#  list("EGAD00001002287", list(
#    list("EGAR00001270598_141030_SN935_0211_B_C5R24ACXX_ACTGAT_3_read1.fastq.gz", "EGAR00001270598_141030_SN935_0211_B_C5R24ACXX_ACTGAT_3_read2.fastq.gz", "EGAR00001270598_141030_SN935_0211_B_C5R24ACXX_ACTGAT_3" )
#  )
#  ),
  list("EGAD00001002315", list(
    list("EGAR00001356286_150616_SN546_0278_B_C794GACXX_ACTTGA_2_read1.fastq.gz", "EGAR00001356286_150616_SN546_0278_B_C794GACXX_ACTTGA_2_read2.fastq.gz", "EGAR00001356286_150616_SN546_0278_B_C794GACXX_ACTTGA_2", 42),
    list("EGAR00001356280_150616_SN546_0278_B_C794GACXX_AGTTCC_5_read1.fastq.gz", "EGAR00001356280_150616_SN546_0278_B_C794GACXX_AGTTCC_5_read2.fastq.gz", "EGAR00001356280_150616_SN546_0278_B_C794GACXX_AGTTCC_5", 43),
    list("EGAR00001384199_151008_SN935_0021_A_C7ANWACXX_TGACCA_5_read1.fastq.gz", "EGAR00001384199_151008_SN935_0021_A_C7ANWACXX_TGACCA_5_read1.fastq.gz", "EGAR00001384199_151008_SN935_0021_A_C7ANWACXX_TGACCA_5", 44),
    list("EGAR00001356283_150616_SN546_0278_B_C794GACXX_GGCTAC_3_read1.fastq.gz", "EGAR00001356283_150616_SN546_0278_B_C794GACXX_GGCTAC_3_read2.fastq.gz", "EGAR00001356283_150616_SN546_0278_B_C794GACXX_GGCTAC_3", 45),
    list("EGAR00001374871_150922_SN935_0020_B_C788LACXX_AGTCAA_8_read1.fastq.gz", "EGAR00001374871_150922_SN935_0020_B_C788LACXX_AGTCAA_8_read2.fastq.gz", "EGAR00001374871_150922_SN935_0020_B_C788LACXX_AGTCAA_8", 46),
    list("EGAR00001356287_150616_SN546_0278_B_C794GACXX_TTAGGC_1_read1.fastq.gz", "EGAR00001356287_150616_SN546_0278_B_C794GACXX_TTAGGC_1_read2.fastq.gz", "EGAR00001356287_150616_SN546_0278_B_C794GACXX_TTAGGC_1", 47)
  )
  ),
  list("EGAD00001002320", list(
    list("EGAR00001384202_151008_SN935_0021_A_C7ANWACXX_AGTCAA_3_read1.fastq.gz", "EGAR00001384202_151008_SN935_0021_A_C7ANWACXX_AGTCAA_3_read2.fastq.gz", "EGAR00001384202_151008_SN935_0021_A_C7ANWACXX_AGTCAA_3", 48)
  )
  ),
#  list("EGAD00001002321", list(
#    list("EGAR00001074529_s_120913_7.CGTACG.read1.fastq.gz", "EGAR00001074529_s_120913_7.CGTACG.read2.fastq.gz", "EGAR00001074529_s_120913_7.CGTACG" ),
#    list("EGAR00001074532_s_120913_6.CAGATC.read1.fastq.gz", "EGAR00001074532_s_120913_6.CAGATC.read2.fastq.gz", "EGAR00001074532_s_120913_6.CAGATC" ),
#    list("EGAR00001094080_130214_SN935_0147_B_D1T9CACXX_GAGTGG_7_read1.fastq.gz", "EGAR00001094080_130214_SN935_0147_B_D1T9CACXX_GAGTGG_7_read2.fastq.gz", "EGAR00001094080_130214_SN935_0147_B_D1T9CACXX_GAGTGG_7" ),
#    list("EGAR00001384203_151008_SN935_0021_A_C7ANWACXX_ACAGTG_3_read1.fastq.gz", "EGAR00001384203_151008_SN935_0021_A_C7ANWACXX_ACAGTG_3_read2.fastq.gz", "EGAR00001384203_151008_SN935_0021_A_C7ANWACXX_ACAGTG_3" )
#  )
#  ),
  list("EGAD00001002347", list(
    list("EGAR00001074528_s_120913_7.GTCCGC.read1.fastq.gz", "EGAR00001074528_s_120913_7.GTCCGC.read2.fastq.gz", "EGAR00001074528_s_120913_7.GTCCGC", 49)
  )
  ),
  list("EGAD00001002349", list(
    list("EGAR00001074534_s_120913_5.GTGGCC.read1.fastq.gz", "EGAR00001074534_s_120913_5.GTGGCC.read2.fastq.gz", "EGAR00001074534_s_120913_5.GTGGCC", 50),
    list("EGAR00001324685_150415_SN546_0275_A_H9A87ADXX_GTGAAA_2_read1.fastq.gz", "EGAR00001324685_150415_SN546_0275_A_H9A87ADXX_GTGAAA_2_read2.fastq.gz", "EGAR00001324685_150415_SN546_0275_A_H9A87ADXX_GTGAAA_2", 51)
  )
  ),
  list("EGAD00001002351", list(
    list("EGAR00001074537_s_120913_4.TAGCTT.read1.fastq.gz", "EGAR00001074537_s_120913_4.TAGCTT.read2.fastq.gz", "EGAR00001074537_s_120913_4.TAGCTT", 52)
  )
  ),
  list("EGAD00001002414", list(
    list("EGAR00001356288_150616_SN546_0278_B_C794GACXX_GATCAG_1_read1.fastq.gz", "EGAR00001356288_150616_SN546_0278_B_C794GACXX_GATCAG_1_read2.fastq.gz", "EGAR00001356288_150616_SN546_0278_B_C794GACXX_GATCAG_1", 53)
  )
  ),
#  list("EGAD00001002426", list(
#    list("EGAR00001356255_150603_SN935_0017_A_C73GHACXX_GTGAAA_7_read1.fastq.gz", "EGAR00001356255_150603_SN935_0017_A_C73GHACXX_GTGAAA_7_read2.fastq.gz", "EGAR00001356255_150603_SN935_0017_A_C73GHACXX_GTGAAA_7" ),
#    list("EGAR00001356275_150616_SN546_0278_B_C794GACXX_GTCCGC_7_read1.fastq.gz", "EGAR00001356275_150616_SN546_0278_B_C794GACXX_GTCCGC_7_read2.fastq.gz", "EGAR00001356275_150616_SN546_0278_B_C794GACXX_GTCCGC_7" ),
#    list("EGAR00001356276_150616_SN546_0278_B_C794GACXX_ACTGAT_7_read1.fastq.gz", "EGAR00001356276_150616_SN546_0278_B_C794GACXX_ACTGAT_7_read2.fastq.gz", "EGAR00001356276_150616_SN546_0278_B_C794GACXX_ACTGAT_7" )
#  )
#  ),
  list("EGAD00001002452", list(
    list("EGAR00001356278_150616_SN546_0278_B_C794GACXX_CTTGTA_6_read1.fastq.gz", "EGAR00001356278_150616_SN546_0278_B_C794GACXX_CTTGTA_6_read2.fastq.gz", "EGAR00001356278_150616_SN546_0278_B_C794GACXX_CTTGTA_6", 54),
    list("EGAR00001356284_150616_SN546_0278_B_C794GACXX_CAGATC_3_read1.fastq.gz", "EGAR00001356284_150616_SN546_0278_B_C794GACXX_CAGATC_3_read2.fastq.gz", "EGAR00001356284_150616_SN546_0278_B_C794GACXX_CAGATC_3", 55),
    list("EGAR00001356277_150616_SN546_0278_B_C794GACXX_GCCAAT_6_read1.fastq.gz", "EGAR00001356277_150616_SN546_0278_B_C794GACXX_GCCAAT_6_read2.fastq.gz", "EGAR00001356277_150616_SN546_0278_B_C794GACXX_GCCAAT_6", 56)
  )
  ),
  list("EGAD00001002456", list(
    list("EGAR00001074538_s_120913_4.CGATGT.read1.fastq.gz", "EGAR00001074538_s_120913_4.CGATGT.read2.fastq.gz", "EGAR00001074538_s_120913_4.CGATGT", 57),
    list("EGAR00001356258_150603_SN935_0017_A_C73GHACXX_ACAGTG_5_read1.fastq.gz", "EGAR00001356258_150603_SN935_0017_A_C73GHACXX_ACAGTG_5_read2.fastq.gz", "EGAR00001356258_150603_SN935_0017_A_C73GHACXX_ACAGTG_5", 58)
  )
  ),
  list("EGAD00001002469", list(
    list("EGAR00001324688_150415_SN546_0275_A_H9A87ADXX_ACTGAT_1_read1.fastq.gz", "EGAR00001324688_150415_SN546_0275_A_H9A87ADXX_ACTGAT_1_read2.fastq.gz", "EGAR00001324688_150415_SN546_0275_A_H9A87ADXX_ACTGAT_1", 59),
    list("EGAR00001074536_s_120913_5.CTTGTA.read1.fastq.gz", "EGAR00001074536_s_120913_5.CTTGTA.read2.fastq.gz", "EGAR00001074536_s_120913_5.CTTGTA", 60)
  )
  ),
  list("EGAD00001002476", list(
    list("EGAR00001408395_151008_SN935_0021_A_C7ANWACXX_TGACCA_2_read1.fastq.gz", "EGAR00001408395_151008_SN935_0021_A_C7ANWACXX_TGACCA_2_read2.fastq.gz", "EGAR00001408395_151008_SN935_0021_A_C7ANWACXX_TGACCA_2", 61),
    list("EGAR00001074530_s_120913_7.ACTGAT.read1.fastq.gz", "EGAR00001074530_s_120913_7.ACTGAT.read2.fastq.gz", "EGAR00001074530_s_120913_7.ACTGAT", 62),
    list("EGAR00001356281_150616_SN546_0278_B_C794GACXX_GTGAAA_4_read1.fastq.gz", "EGAR00001356281_150616_SN546_0278_B_C794GACXX_GTGAAA_4_read2.fastq.gz", "EGAR00001356281_150616_SN546_0278_B_C794GACXX_GTGAAA_4", 63)
  )
  ),
  list("EGAD00001002482", list(
    list("EGAR00001074539_s_120913_4.AGTCAA.read1.fastq.gz", "EGAR00001074539_s_120913_4.AGTCAA.read2.fastq.gz", "EGAR00001074539_s_120913_4.AGTCAA", 64)
  )
  )
)	

dataFolderBP = "E:/BulkProfiles/Blueprint/FASTQ/EGA_download_client"


for (proj in 1:length(BPFileInfo)) {
#for (proj in 1) {
  projName = BPFileInfo[[proj]][[1]]
  samples = BPFileInfo[[proj]][[2]]
  for (sampInd in 1:length(samples)) {
    sampleFiles = samples[[sampInd]];
    outputFolder = paste0(dataFolderBP, "/", projName, "/", sampleFiles[[3]])
    
    ## Read in abundance file from kallisto
    abundance = read.table(file=paste0(outputFolder, "/abundance.tsv"), header=T)
    head(abundance)
    
    colnames(abundance) = c("ENST", "length", "eff_length", "est_counts", "tpm")
    head(abundance)
    head(gene_names)
    
    # Align
    aligned = merge(abundance, gene_names, by = "ENST")
    head(aligned)
    
    # Get rid of everything but est_counts and tpm
    aligned = aligned[ ,c("ENST", "ENSG", "HGNC", "est_counts", "tpm")]
    head(aligned)
    
    # Aggregate
    agg = aggregate(. ~ HGNC, data = aligned, sum)
    agg = agg[,c('HGNC','est_counts','tpm')]
    head(agg)
    
    # Write agg to file
    write.table(agg, file=paste0(outputFolder, "/abundance_hgnc.txt"),
                row.names=F, sep="\t")
  }
}

# Make a big data frame with all (ENCODE, FANTOM5, BLUEPRINT, GSE 51984) samples. HGNC gene as rows and subject/sample as column

# Read in the first sample
count_df = read.table(file=paste0(dataFolderEF, "/output_1/abundance_hgnc.txt"), header=T, sep="\t")
# Set hgnc as rownames
rownms = count_df[ ,1]
row.names(count_df) = rownms
head(count_df)
count_df = count_df[ ,-1]
head(count_df)
dim(count_df)

tpm_df = count_df[, "tpm", drop=FALSE]
count_df = count_df[, "est_counts", drop=FALSE]

head(count_df)

# Read in the rest
for ( i in 2:25) {
  d = read.table(file=paste0(dataFolderEF, "/output_", i, "/abundance_hgnc.txt"), header=T, sep="\t");
  count_df = cbind(count_df, d[ ,"est_counts"] )
  tpm_df = cbind(tpm_df, d[ ,"tpm"] )
}

colnames(count_df) = seq(1,25,1)
head(count_df)
colnames(tpm_df) = seq(1,25,1)
head(tpm_df)

head(count_df)
dim(count_df)
head(tpm_df)
dim(tpm_df)


#add the blueprint data as well
for (proj in 1:length(BPFileInfo)) {
  #for (proj in 1) {
  projName = BPFileInfo[[proj]][[1]]
  samples = BPFileInfo[[proj]][[2]]
  for (sampInd in 1:length(samples)) {
    sampleFiles = samples[[sampInd]];
    outputFolder = paste0(dataFolderBP, "/", projName, "/", sampleFiles[[3]])
    
    d = read.table(file=paste0(outputFolder, "/abundance_hgnc.txt"), header=T, sep="\t");
    count_df = cbind(count_df, d[ ,"est_counts"] )
    tpm_df = cbind(tpm_df, d[ ,"tpm"] )
    colnames(count_df)[[length(colnames(count_df))]] = sampleFiles[[4]];
    colnames(tpm_df)[[length(colnames(tpm_df))]] = sampleFiles[[4]];
  }
}

#now add the GSE51984 samples. These are not processed the same way, so the genes
#need to be syncronized

dataFolderGSE51984 = "C:/Work/MatlabCode/components/SCLib/ImportableData/GSE51984"

fn = list("GSM1256812_B_01_RPKM.txt",
          "GSM1256813_B_02_RPKM.txt",
          "GSM1256814_B_03_RPKM.txt",
          "GSM1256815_B_04_RPKM.txt",
          "GSM1256816_B_05_RPKM.txt",
          "GSM1256828_T_01_RPKM.txt",
          "GSM1256829_T_02_RPKM.txt",
          "GSM1256830_T_03_RPKM.txt",
          "GSM1256831_T_04_RPKM.txt",
          "GSM1256832_T_05_RPKM.txt"
  )

count_df2 = read.table(file=paste0(dataFolderGSE51984, "/", fn[[1]]), header=F, sep="\t")
colnames(count_df2)[[4]] = "HGNC";

tpm_df2 = count_df2[, c(4,5), drop=FALSE] #so, TPM is FPKM, this will be fixed later
count_df2 = count_df2[, c(4,6), drop=FALSE]

lastIndex = 64

colnames(count_df2)[[length(colnames(count_df2))]] = lastIndex+1;
colnames(tpm_df2)[[length(colnames(tpm_df2))]] = lastIndex+1;

head(count_df2)


for (i in 2:10) {
  d = read.table(file=paste0(dataFolderGSE51984, "/", fn[[i]]), header=F, sep="\t")
  count_df2 = cbind(count_df2, d[, 6] )
  tpm_df2 = cbind(tpm_df2, d[ , 5] )
  colnames(count_df2)[[length(colnames(count_df2))]] = lastIndex+i;
  colnames(tpm_df2)[[length(colnames(tpm_df2))]] = lastIndex+i;
}

#some genes exist at two places in the genome and thus have two rows; aggregate those
# Aggregate
count_df2 = aggregate(. ~ HGNC, data = count_df2, sum)
tpm_df2 = aggregate(. ~ HGNC, data = tpm_df2, sum)
#head(aggCount)
#dim(count_df2)
#dim(aggCount)


# Set hgnc as rownames
rownms = count_df2[[colnames(count_df2)[1]]]
row.names(count_df2) = rownms
count_df2 = count_df2[ ,-1]
head(count_df2)
dim(count_df2)
row.names(tpm_df2) = tpm_df2[[colnames(tpm_df2)[1]]]
tpm_df2 = tpm_df2[ ,-1]
head(tpm_df2)

library("data.table")

#must use data.table; data.frame, matrix and similar are super slow and takes up too much memory,
#which is super weird, this shouldn't be a problem at all!
c1 = setDT(count_df, keep.rownames=TRUE)
c2 = setDT(count_df2, keep.rownames=TRUE)

#now merge the datasets, throwing away all genes that does not exist in both datasets
mergedCounts = merge(c1, c2)
dim(mergedCounts)

counts = as.data.frame(mergedCounts)
row.names(counts) = counts[, 1];
counts = counts[,-1]
counts = data.matrix(counts)

c1 = setDT(tpm_df, keep.rownames=TRUE)
c2 = setDT(tpm_df2, keep.rownames=TRUE)

#now merge the datasets, throwing away all genes that does not exist in both datasets
#mergedCounts = merge(count_df, count_df2) - takes forever and then runs out of memory
mergedTPMs = merge(c1, c2)
dim(mergedTPMs)

tpms = as.data.frame(mergedTPMs)
row.names(tpms) = tpms[, 1];
tpms = tpms[,-1]

tpms = data.matrix(tpms)



#rescale to TPM
tpm <- function(data){
  for(i in 1:dim(data)[2]){
    data[,i] <- data[,i] * 1e6 / sum(data[,i])
  }
  data[is.na(data)] <- 0
  return(data)
}

tpms = tpm(tpms)
#apply(tpms, 2, sum)

#now save to disk:
write.table(counts, file=paste0(dataFolderEF, "/countsMatrix.txt"),
            row.names=T, sep="\t")

write.table(tpms, file=paste0(dataFolderEF, "/tpmMatrix.txt"),
            row.names=T, sep="\t")

#also save a TMM normalized dataset
#library("edgeR", lib.loc="~/R/win-library/3.5")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("limma")
#BiocManager::install("edgeR")


#data = read.csv(inFilename, header=TRUE, row.names=1, sep = "\t", quote = "\"'", check.names=FALSE);
#using TMM from edgeR:
normFactors <-  calcNormFactors(tpms)
#for some weird reason I need to transpose the data matrix back and forth to get the
#row wise multiplication to work...
A = t(tpms);
B = A / normFactors;
  
tmms = t(B);
  
write.table(tmms, file=paste0(dataFolderEF, "/tmmMatrix.txt"),
            row.names=T, sep="\t")

#install.packages("R.matlab")
#library("R.matlab")

writeMat(paste0(dataFolderEF, "/tmms.mat"), tmmMat = tmms, tmmGenes = row.names(tmms), tmmSampIds = colnames(tmms))
writeMat(paste0(dataFolderEF, "/tpms.mat"), tpmMat = tpms, tpmGenes = row.names(tpms), tpmSampIds = colnames(tpms))
writeMat(paste0(dataFolderEF, "/counts.mat"), countsMat = counts, countsGenes = row.names(counts), tpmSampIds = colnames(counts))

#experiment with quantile normalization as well:
#BiocManager::install("preprocessCore")
library("preprocessCore")

quantileNorms = normalize.quantiles(tpms)
colSums(quantileNorms)
row.names(quantileNorms) = row.names(tpms);
colnames(quantileNorms) = colnames(tpms)
head(quantileNorms)
writeMat(paste0(dataFolderEF, "/qns.mat"), qnsMat = quantileNorms, qnsGenes = row.names(quantileNorms), qnsSampIds = colnames(quantileNorms))

#also try RUV:
#BiocManager::install("RUVSeq")
#library("RUVSeq")

#read the design matrix
#install.packages("xlsx")
library("xlsx")
dm <- read.xlsx("C:/Work/R/RNASeqCTProfileEval/DesignMatrix.xlsx", sheetName = "DesignMatrix")
cellTypes = as.numeric(dm[9, 2:75])
labs = as.numeric(dm[2, 2:75])
subCellTypes = as.numeric(dm[7, 2:75])
tissues = as.numeric(dm[3, 2:75])


#get housekeeping genes
pathHkGenes = "C:/Work/MatlabCode/components/SCLib/ImportableData/HK_genes.txt"
hkGenes = read.table(pathHkGenes)
hkGenes = hkGenes[,1]


library(RUVSeq)
filter <- apply(counts, 1, function(x) length(x[x>5])>=2)
filtered <- counts[filter,]

#genes <- rownames(filtered)[grep("^ENS", rownames(filtered))]
#spikes <- rownames(filtered)[grep("^ERCC", rownames(filtered))]

x <- as.factor(cellTypes)

library(RColorBrewer)
colors <- brewer.pal(3, "Set2")
plotRLE(counts, outline=FALSE, ylim=c(-4, 4), col=as.numeric(cellTypes))
plotPCA(counts, col=as.numeric(cellTypes), cex=1.2)

uc = betweenLaneNormalization(filtered, which="upper")
plotRLE(uc, outline=FALSE, ylim=c(-4, 4), col=as.numeric(cellTypes))
plotPCA(uc, col=as.numeric(cellTypes), cex=1.2)

plotRLE(tmms, outline=FALSE, ylim=c(-4, 4), col=as.numeric(cellTypes))
plotPCA(tmms, col=as.numeric(cellTypes), cex=1.2)


plotRLE(tpms, outline=FALSE, ylim=c(-4, 4), col=as.numeric(cellTypes))


#it seems that upper quartile normalization does not do a very good job, compared to TMM?


RUVInput = newSeqExpressionSet(as.matrix(uc),
                           phenoData = data.frame(x, row.names=colnames(uc)))

#RUVInput = newSeqExpressionSet(as.matrix(uc),
#                               phenoData = data.frame(x, row.names=colnames(tmms)))


hkGenes2 = intersect(hkGenes, row.names(uc))

RUVed = RUVg(RUVInput, hkGenes2, k=1)
plotPCA(RUVed, col=cellTypes, cex=1.2)
plotPCA(RUVed, col=labs, cex=1.2)

#test combat

#BiocManager::install("sva")
#library("sva")
#library(bladderbatch)
#data(bladderdata)
#library(pamr)
#library(limma)

#first, filter lowly expressed genes (based on counts)
filter <- apply(counts, 1, function(x) length(x[x>5])>=2)
filttmms = tmms[filter,]

#create log-transformed data from TMM normalization
#logtmms = unlist(as.data.frame(log2(filttmms + 1)))
logtmms = as.matrix(log2(filttmms + 1))
ctVar = cellTypes - 1;#0 or 1 depending on b cell or t cell
batch = labs;

modcombat = model.matrix(~1 + ctVar, data=as.data.frame(ctVar))

#modcombat2 = model.matrix(~1, data=pheno)
#edata = exprs(bladderEset)

#batch2 = pheno$batch
#combat_edata2 = ComBat(dat=edata, batch=batch2, mod=modcombat2, par.prior=TRUE, prior.plots=FALSE)


combat_edata = ComBat(dat=logtmms, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

#try pca on corrected data
plotPCA((2^combat_edata)-1, col=as.numeric(cellTypes), cex=1.2)
plotPCA((2^combat_edata)-1, col=as.numeric(subCellTypes), cex=1.2)
plotPCA((2^combat_edata)-1, col=as.numeric(tissues), cex=1.2)
plotPCA((2^combat_edata)-1, col=as.numeric(labs), cex=1.2)

plotPCA(tmms, col=as.numeric(labs), cex=1.2)

