%This file uses SingleCellToolbox to load public single-cell datasets and creates 
%pooled gene expression profiles from B and T cells from them. The profiles are exported
%to a text file for importing into R. The genes are synchronized between the datasets, meaning
%that all genes that don't exist in all datasets are lost.


%% First synchronize the genes

[lct, lch] = DsLC.get();
[~,~,gse112845cd8] = DsGSE112845.get();
dsList = {DsHcaCB.get(), lct, lch, DsPbmc68k.get(), DsTCD4Mem.get(), DsB10k.get(), gse112845cd8, DsMel.get()};

dsList = SynchronizeGenes(dsList, [], true);

hca = dsList{1};
lct = dsList{2};
lch = dsList{3};
pbmc68k = dsList{4};
tcd4mem = dsList{5};
b10k = dsList{6};
gse112845cd8 = dsList{7};
mel = dsList{8};


%% HCA TCells and BCells from 8 patients, in total 16 samples

hcat = hca.cellSubset(hca.cellType == Celltype.TCellCD4Pos | hca.cellType == Celltype.TCellCD8Pos | hca.cellType == Celltype.TCellReg |  hca.cellType == Celltype.TCell );
hcab = hca.cellSubset(hca.cellType == Celltype.BCell);

hcat1 = hcat.cellSubset(strcmp(hcat.sampleIds,'CB1'));
hcat2 = hcat.cellSubset(strcmp(hcat.sampleIds,'CB2'));
hcat3 = hcat.cellSubset(strcmp(hcat.sampleIds,'CB3'));
hcat4 = hcat.cellSubset(strcmp(hcat.sampleIds,'CB4'));
hcat5 = hcat.cellSubset(strcmp(hcat.sampleIds,'CB5'));
hcat6 = hcat.cellSubset(strcmp(hcat.sampleIds,'CB6'));
hcat7 = hcat.cellSubset(strcmp(hcat.sampleIds,'CB7'));
hcat8 = hcat.cellSubset(strcmp(hcat.sampleIds,'CB8'));

hcab1 = hcab.cellSubset(strcmp(hcab.sampleIds,'CB1'));
hcab2 = hcab.cellSubset(strcmp(hcab.sampleIds,'CB2'));
hcab3 = hcab.cellSubset(strcmp(hcab.sampleIds,'CB3'));
hcab4 = hcab.cellSubset(strcmp(hcab.sampleIds,'CB4'));
hcab5 = hcab.cellSubset(strcmp(hcab.sampleIds,'CB5'));
hcab6 = hcab.cellSubset(strcmp(hcab.sampleIds,'CB6'));
hcab7 = hcab.cellSubset(strcmp(hcab.sampleIds,'CB7'));
hcab8 = hcab.cellSubset(strcmp(hcab.sampleIds,'CB8'));

%figure
%histogram(categorical(hca.sampleIds))

%% LC : use pat 3, 4 and 5 from the cancer and all patients for the healthy
%tissue; there are fewer cells there

lctpat3 = lct.cellSubset(strcmp(lct.sampleIds,'3'));
lctpat4 = lct.cellSubset(strcmp(lct.sampleIds,'4'));
lctpat5 = lct.cellSubset(strcmp(lct.sampleIds,'5'));

lctpat3t = lctpat3.cellSubset(lctpat3.cellType == Celltype.TCellCD4Pos | lctpat3.cellType == Celltype.TCellCD8Pos | lctpat3.cellType == Celltype.TCellReg |  lctpat3.cellType == Celltype.TCell );
lctpat4t = lctpat4.cellSubset(lctpat4.cellType == Celltype.TCellCD4Pos | lctpat4.cellType == Celltype.TCellCD8Pos | lctpat4.cellType == Celltype.TCellReg |  lctpat4.cellType == Celltype.TCell );
lctpat5t = lctpat5.cellSubset(lctpat5.cellType == Celltype.TCellCD4Pos | lctpat5.cellType == Celltype.TCellCD8Pos | lctpat5.cellType == Celltype.TCellReg |  lctpat5.cellType == Celltype.TCell );

lctpat3b = lctpat3.cellSubset(lctpat3.cellType == Celltype.BCell);
lctpat4b = lctpat4.cellSubset(lctpat4.cellType == Celltype.BCell);
lctpat5b = lctpat5.cellSubset(lctpat5.cellType == Celltype.BCell);

lcht = lch.cellSubset(lch.cellType == Celltype.TCellCD4Pos | lch.cellType == Celltype.TCellCD8Pos | lch.cellType == Celltype.TCellReg |  lch.cellType == Celltype.TCell );
lchb = lch.cellSubset(lch.cellType == Celltype.BCell);

%% PBMC68k : one sample of each, this is just one patient
pbmc68kt = pbmc68k.cellSubset(pbmc68k.cellType == Celltype.TCellCD4Pos | pbmc68k.cellType == Celltype.TCellCD8Pos | pbmc68k.cellType == Celltype.TCellReg |  pbmc68k.cellType == Celltype.TCell );
pbmc68kb = pbmc68k.cellSubset(pbmc68k.cellType == Celltype.BCell);

%% TCD4Mem
%nothing needs to be done


%% B10k
%nothing needs to be done



%% GSE112845
%nothing needs to be done

%% Mel
melt = mel.cellSubset(mel.cellType == Celltype.TCellCD4Pos | mel.cellType == Celltype.TCellCD8Pos | mel.cellType == Celltype.TCellReg |  mel.cellType == Celltype.TCell );
melb = mel.cellSubset(mel.cellType == Celltype.BCell);


%% Export them all to a Samples object
dss = { hcat1, hcat2, hcat3, hcat4, hcat5, hcat6, hcat7, hcat8, hcab1, hcab2, hcab3, hcab4, hcab5, hcab6, hcab7, hcab8, ...
        lctpat3t, lctpat4t, lctpat5t, lctpat3b, lctpat4b, lctpat5b, lcht, lchb, pbmc68kt, pbmc68kb, tcd4mem, b10k, gse112845cd8, melt, melb};
numSets = size(dss,2);
numGenes = size(hcat1.data,1);

samp = Samples;
samp.data = zeros(numGenes,numSets);
samp.sampleIds = {'hcat1', 'hcat2', 'hcat3', 'hcat4', 'hcat5', 'hcat6', 'hcat7', 'hcat8', 'hcab1', 'hcab2', 'hcab3', 'hcab4', 'hcab5', ...
                  'hcab6', 'hcab7', 'hcab8', 'lctpat3t', 'lctpat4t', 'lctpat5t', 'lctpat3b', 'lctpat4b', 'lctpat5b', 'lcht', 'lchb', ...
                  'pbmc68kt', 'pbmc68kb', 'tcd4mem', 'b10k', 'gse112845cd8', 'melt', 'melb'};
samp.genes = hcat1.genes;

for i = 1:numSets
    ds = dss{1,i};
    %We do not TPM normalize here, since we want to keep the counts information
	%to be able to estimate the dispersion in TMM. Normalization is done in R later.
    datasum = sum(ds.data,2);
    samp.data(:,i) = datasum;
    %also count the number of cells per sample
end

samp.writeToTextFile('data/scProfiles.txt');

%Now, also create single-cell profiles for cibersort. For practical
%reasons (due to storage limitations in cibersort), we do not produce more 
%than 100 cells per cell type and dataset. Since this is far too few to get a stable expression, we
%create "cells" as pools of many cells

GenCibersortProfiles(mel, 'data/deconv/sc/melProfiles', 100, 100);
GenCibersortProfiles(mel, 'data/deconv/sc/melProfilesUneven', 20, 80);
GenCibersortProfiles(hca, 'data/deconv/sc/hcaProfiles', 100, 100);
GenCibersortProfiles(lct, 'data/deconv/sc/lctProfiles', 100, 100);
GenCibersortProfiles(lch, 'data/deconv/sc/lchProfiles', 100, 100);
GenCibersortProfiles(pbmc68k, 'data/deconv/sc/pbmc68kProfiles', 100, 100);

%get how many cells that were used for each profile

ds = mel;
sum(ds.cellType == Celltype.BCell)
%512
sum(ds.cellType == Celltype.TCellCD4Pos | ds.cellType == Celltype.TCellCD8Pos | ds.cellType == Celltype.TCellReg |  ds.cellType == Celltype.TCell)
%2040

ds = hca;
sum(ds.cellType == Celltype.BCell)
%35910
sum(ds.cellType == Celltype.TCellCD4Pos | ds.cellType == Celltype.TCellCD8Pos | ds.cellType == Celltype.TCellReg |  ds.cellType == Celltype.TCell)
%167612

ds = lct;
sum(ds.cellType == Celltype.BCell)
%4509
sum(ds.cellType == Celltype.TCellCD4Pos | ds.cellType == Celltype.TCellCD8Pos | ds.cellType == Celltype.TCellReg |  ds.cellType == Celltype.TCell)
%18306

ds = lch;
sum(ds.cellType == Celltype.BCell)
%297
sum(ds.cellType == Celltype.TCellCD4Pos | ds.cellType == Celltype.TCellCD8Pos | ds.cellType == Celltype.TCellReg |  ds.cellType == Celltype.TCell)
%4864

ds = pbmc68k;
sum(ds.cellType == Celltype.BCell)
%5908
sum(ds.cellType == Celltype.TCellCD4Pos | ds.cellType == Celltype.TCellCD8Pos | ds.cellType == Celltype.TCellReg |  ds.cellType == Celltype.TCell)
%48657

%get the number of cells in each pooled sample (entered into the design matrix excel sheet)
length(hcat1.cellIds)
length(hcat2.cellIds)
length(hcat3.cellIds)
length(hcat4.cellIds)
length(hcat5.cellIds)
length(hcat6.cellIds)
length(hcat7.cellIds)
length(hcat8.cellIds)

length(hcab1.cellIds)
length(hcab2.cellIds)
length(hcab3.cellIds)
length(hcab4.cellIds)
length(hcab5.cellIds)
length(hcab6.cellIds)
length(hcab7.cellIds)
length(hcab8.cellIds)

length(lctpat3t.cellIds)
length(lctpat4t.cellIds)
length(lctpat5t.cellIds)
length(lctpat3b.cellIds)
length(lctpat4b.cellIds)
length(lctpat5b.cellIds)
length(lcht.cellIds)
length(lchb.cellIds)

length(pbmc68kt.cellIds)
length(pbmc68kb.cellIds)
length(tcd4mem.cellIds)
length(b10k.cellIds)
length(gse112845cd8.cellIds)
length(melt.cellIds)
length(melb.cellIds)


