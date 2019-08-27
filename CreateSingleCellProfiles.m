%% First synchronize the genes

[lct, lch] = DsLC.get();
[~,~,gse112845cd8] = DsGSE112845.get();
dsList = {DsHcaCB.get(), lct, lch, DsPbmc68k.get(), DsTCD4Mem.get(), DsB10k.get(), gse112845cd8};

dsList = SynchronizeGenes(dsList, [], true);

hca = dsList{1};
lct = dsList{2};
lch = dsList{3};
pbmc68k = dsList{4};
tcd4mem = dsList{5};
b10k = dsList{6};
gse112845cd8 = dsList{7};



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


%% Export them all to a Samples object
%Synchronize the genes:
dss = { hcat1, hcat2, hcat3, hcat4, hcat5, hcat6, hcat7, hcat8, hcab1, hcab2, hcab3, hcab4, hcab5, hcab6, hcab7, hcab8, ...
        lctpat3t, lctpat4t, lctpat5t, lctpat3b, lctpat4b, lctpat5b, lcht, lchb, pbmc68kt, pbmc68kb, tcd4mem, b10k, gse112845cd8};
numSets = size(dss,2);
numGenes = size(hcat1.data,1);

samp = Samples;
samp.data = zeros(numGenes,numSets);
samp.sampleIds = {'hcat1', 'hcat2', 'hcat3', 'hcat4', 'hcat5', 'hcat6', 'hcat7', 'hcat8', 'hcab1', 'hcab2', 'hcab3', 'hcab4', 'hcab5', ...
                  'hcab6', 'hcab7', 'hcab8', 'lctpat3t', 'lctpat4t', 'lctpat5t', 'lctpat3b', 'lctpat4b', 'lctpat5b', 'lcht', 'lchb', ...
                  'pbmc68kt', 'pbmc68kb', 'tcd4mem', 'b10k', 'gse112845cd8'};
samp.genes = hcat1.genes;

totcounts = zeros(1,numSets);

for i = 1:numSets
    ds = dss{1,i};
    %we do not tpm normalize until after summing up all cells, since we
    %want the counts to truly represent the noise. This should not matter
    %much
    datasum = sum(ds.data,2);
    samp.data(:,i) = TPM(datasum);
    totcounts(1,i) = sum(datasum,1);
end

samp.writeToTextFile('scProfiles.txt');
dlmwrite('scTotCounts.txt',totcounts,'\t');


