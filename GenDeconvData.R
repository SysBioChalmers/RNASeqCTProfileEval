

dataFolder = "C:/Work/R/RNASeqCTProfileEval/"
source(paste0(dataFolder, "FigureHelpFunc.R"))



#############################
# load data
#############################
library(dplyr)
library(tibble)

tpmScAndBulk = readRDS(paste0(dataFolder, "data/tpmScAndBulk.RDS"))
tmmScAndBulk = readRDS(paste0(dataFolder, "data/tmmScAndBulk.RDS"))
quantileScAndBulk = readRDS(paste0(dataFolder, "data/quantileScAndBulk.RDS"))
bcScAndBulk = readRDS(paste0(dataFolder, "data/bcScAndBulk.RDS"))

tpmScAndBulkNonFilt = readRDS(paste0(dataFolder, "data/tpmScAndBulkNonFilt.RDS"))
tmmScAndBulkNonFilt = readRDS(paste0(dataFolder, "data/tmmScAndBulkNonFilt.RDS"))
bcScAndBulk = readRDS(paste0(dataFolder, "data/quantileScAndBulkNonFilt.RDS"))

cellTypes = readRDS(paste0(dataFolder, "data/cellTypes.RDS"))
labs = readRDS(paste0(dataFolder, "data/labs.RDS"))
subCellTypes = readRDS(paste0(dataFolder, "data/subCellTypes.RDS"))
tissues = readRDS(paste0(dataFolder, "data/tissues.RDS"))
individual = readRDS(paste0(dataFolder, "data/individual.RDS"))
techRepl = readRDS(paste0(dataFolder, "data/techRepl.RDS"))
scOrBulk = readRDS(paste0(dataFolder, "data/scOrBulk.RDS"))


profilesAll = tmmScAndBulkNonFilt
cellTypeNames = rep("TCells", length(cellTypes))
cellTypeNames[cellTypes == 1] = "BCells"
#colnames(profilesAll) = cellTypeNames


#create bulk cell type profiles from labs (skip lab 3, has no t cells)
##########################################
CreateCellTypeProfiles <- function(sel, profileId, exprData) {
  #first, check that we have at least 2 of each type - if not, duplicate the single column
  #to create two, cibersortx fails otherwise
  origInd = which(sel)
  ind = origInd
  ct = cellTypes[sel]
  if (sum(ct == 1) == 1) {
    ind = c(ind, origInd[ct==1])
  }
  if (sum(ct == 2) == 1) {
    ind = c(ind, origInd[ct==2])
  }
  
  #write the profile file
  profilesSpec = exprData[, ind]
  cellTypesSpec = cellTypes[ind]
  cellTypeNamesSpec = cellTypeNames[ind]
  
  profilesSpec = cbind(GeneSymbol = row.names(profilesSpec), profilesSpec)
  
  write.table(profilesSpec, paste0(dataFolder, "data/deconv/", profileId, "profiles.txt"), row.names = F, col.names = c("GeneSymbol",cellTypeNamesSpec), quote = F, sep="\t")
  
  #write the phenotype file
  bcellRow = rep(0, length(cellTypesSpec))
  
  bcellRow[cellTypesSpec == 1] = 1
  bcellRow[cellTypesSpec != 1] = 2
  tcellRow = rep(0, sum(sel))
  tcellRow[cellTypesSpec == 2] = 1
  tcellRow[cellTypesSpec != 2] = 2
  
  phenTypeSpec = rbind(bcellRow, tcellRow)
  row.names(phenTypeSpec) = c("BCells", "TCells")
  write.table(phenTypeSpec, paste0(dataFolder, "data/deconv/", profileId, "phenotype.txt"), col.names = F, row.names = T, quote = F, sep="\t")
}

CreateCellTypeProfiles(labs == 1, "lab1", profilesAll)
CreateCellTypeProfiles(labs == 4, "lab4", profilesAll)
CreateCellTypeProfiles(labs == 5, "lab5", profilesAll)
CreateCellTypeProfiles(labs == 6, "lab6", profilesAll)
CreateCellTypeProfiles(labs == 7, "lab7", profilesAll)
CreateCellTypeProfiles(labs == 1, "lab1TPM", tpmScAndBulkNonFilt)
CreateCellTypeProfiles(labs == 4, "lab4TPM", tpmScAndBulkNonFilt)
CreateCellTypeProfiles(labs == 5, "lab5TPM", tpmScAndBulkNonFilt)
CreateCellTypeProfiles(labs == 6, "lab6TPM", tpmScAndBulkNonFilt)
CreateCellTypeProfiles(labs == 7, "lab7TPM", tpmScAndBulkNonFilt)




#Create bulk mixes:
####################################

mixes = matrix(data = 0, nrow = dim(tmmScAndBulkNonFilt)[1], ncol = 200)#make it big enough, reduce it later
fracB = 0.5
fracT = 1 - fracB
#create mixtures from the rest of the bulk labs, with exactly 50/50 mix (i.e. labs 1, 2, and 5 - lab 3 has no T cells)
#create 34 mixtures for now

#create mix names
mixNames = rep("", 200)
mixNames[1] = "GeneSymbols"
#create mix pairs
mixPairs = matrix(data = 0, nrow = 2, ncol = 200) # B cell sample in row 1, T cell in row 2


sampleIndex = 1

for (lab in c(1,2,4,5)) {

  bind = which(cellTypes == 1 & labs == lab)
  tind = which(cellTypes == 2 & labs == lab)
  
  numB = length(bind)
  numT = length(tind)
  numSamp = max(numB, numT)
  
  #repeat the smaller one until it reaches the length of the longer one
  if (numB < numT) {
    repTimes = ceiling(numT/numB)
    res = rep(bind, repTimes)
    bind = res[1:numSamp]
  } else if (numT < numB) {
    repTimes = ceiling(numB/numT)
    res = rep(tind, repTimes)
    tind = res[1:numSamp]
  }
  #now add the mixed samples
  for (i in 1:numSamp) {
    #so, use tmm data when taking 50% of each, otherwise the real mix could be off. Then 
    #TPM normalize the new sample before exporting
    mixes[, sampleIndex] = tmmScAndBulkNonFilt[, bind[i]] * fracB + tmmScAndBulkNonFilt[, tind[i]] * fracT
    mixes[, sampleIndex] = mixes[, sampleIndex]*(10^6/sum(mixes[, sampleIndex]))
    mixNames[sampleIndex+1] = paste0("Lab", lab, "_B_", bind[i], "_T_", tind[i])
    mixPairs[1,sampleIndex] = bind[i]
    mixPairs[2,sampleIndex] = tind[i]
    sampleIndex = sampleIndex + 1;
  }
}

#reduce the data to what has been written:
mixes = mixes[,1:(sampleIndex-1)]
mixPairs = mixPairs[,1:(sampleIndex-1)]
mixNames = mixNames[1:sampleIndex]

#write mixes:
#add the row names as a separate column to remove the problem with the empty left corner
df = as.data.frame(mixes)
df = cbind(GeneSymbol = row.names(tmmScAndBulkNonFilt), df)
write.table(df, paste0(dataFolder, "data/deconv/bulkMixes.txt"), row.names = F, col.names = mixNames, quote = F, sep="\t")
#write pairs:
write.table(mixPairs, paste0(dataFolder, "data/deconv/bulkMixes_pairs.txt"), row.names = F, col.names = F, quote = F, sep="\t")









#Now create profiles and mixtures for the same lab (lab 4)
################################

#Lab 4 has two individuals - use one for profiles and predict on the others
#Profiles of individual 1
CreateCellTypeProfiles(labs == 4 & individual == 1, "Lab4Internal", tpmScAndBulkNonFilt)

#mixes of individual 3. Create all combinations of B and T cells
mixesInt = matrix(data = 0, nrow = dim(tpmScAndBulkNonFilt)[1], ncol = 9)

mixNamesInt = rep("", 9)
mixNamesInt[1] = "GeneSymbols"
#create mix pairs
mixPairsInt = matrix(data = 0, nrow = 2, ncol = 9) # B cell sample in row 1, T cell in row 2

sampleIndex = 1

bind = which(cellTypes == 1 & individual == 3)
tind = which(cellTypes == 2 & individual == 3)

numB = length(bind)
numT = length(tind)

#now add the mixed samples
for (b in 1:numB) {
  for (t in 1:numT) {
    mixesInt[, sampleIndex] = tmmScAndBulkNonFilt[, bind[b]] * fracB + tmmScAndBulkNonFilt[, tind[t]] * fracT
    mixesInt[, sampleIndex] = mixesInt[, sampleIndex]*(10^6/sum(mixes[, sampleIndex]))
    mixNamesInt[sampleIndex+1] = paste0("Lab", lab, "_B_", bind[b], "_T_", tind[t])
    mixPairsInt[1,sampleIndex] = bind[b]
    mixPairsInt[2,sampleIndex] = tind[t]
    sampleIndex = sampleIndex + 1;
  }
}

#write mixes:
#add the row names as a separate column to remove the problem with the empty left corner
df = as.data.frame(mixesInt)
df = cbind(GeneSymbol = row.names(tpmScAndBulkNonFilt), df)
write.table(df, paste0(dataFolder, "data/deconv/bulkMixesInt.txt"), row.names = F, col.names = mixNamesInt, quote = F, sep="\t")
#write pairs:
write.table(mixPairsInt, paste0(dataFolder, "data/deconv/bulkMixesInt_pairs.txt"), row.names = F, col.names = F, quote = F, sep="\t")




