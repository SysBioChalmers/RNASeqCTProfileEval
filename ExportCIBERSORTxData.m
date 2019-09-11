
%% 1. Create single-cell profiles for CIBERSORTx. Use all in the cancer part of the lc dataset
[lct, ~] = DsLC.get();

unique(CelltypeId2CelltypeName(lct.cellType))
figure
histogram(lct.cellType)

%remove rare cell types for which there will not be eonough cells to pool
lct = lct.cellSubset((lct.cellType ~= Celltype.Dendritic) & (lct.cellType ~= Celltype.Epithelial) & (lct.cellType ~= Celltype.Langerhans) );

%cap to 2000 cells of each type
% Skip this step, it seems cibersort couldn't handle the file. Instead,
% pool a bunch of cells and treat as one cell
%{
cts = unique(lct.cellType);
for i = 1:size(cts,2)
   sel = lct.cellType == cts(1,i);
   num = sum(sel);
   ind = find(sel);
   indind = randsample(num, num-2000);
   allind = 1:size(lct.data,2);
   allind(:,ind(1,indind)) = [];%remove the samples
   lct = lct.cellSubset(allind);
end
%}

%make sure there are no spaces in the names
names = CelltypeId2CelltypeName(lct.cellType);
names(1, lct.cellType == Celltype.BCell) = {'BCell'};
names(1, lct.cellType == Celltype.TCell) = {'TCell'};
names(1, lct.cellType == Celltype.TCellCD4Pos) = {'TCD4'};
names(1, lct.cellType == Celltype.TCellCD8Pos) = {'TCD8'};
names(1, lct.cellType == Celltype.TCellReg) = {'TReg'};
names(1, lct.cellType == Celltype.NKCell) = {'NKCell'};

unique(names)

%now, join 20 cells of the same type and turn into one


toExp = lct;

index = 1;

cts = unique(lct.cellType);
for i = 1:size(cts,2)
    cts(i)
    sel = lct.cellType == cts(1,i);
    num = sum(sel);
    ind = find(sel);
    for j = 1:20:num
        j
        endp = min(j+19,num);
        toExp.data(:,index) = sum(lct.data(:,ind(j:endp)),2);
        toExp.cellType(1,index) = cts(i);
        index = index + 1;
    end
end

toExp = toExp.cellSubset(1:(index-1));


%toExp.cellIds = names; %set cellIds to names to be able to export with cell type as header
%toExp.saveDataTable('CIBERSORTx_scProfiles.txt');

%export via a table doesn't work since the cell types should be heading. So
%we need to write the file manually instead:
tic
[fid,msg] = fopen('CIBERSORTx_scProfiles4.txt','wt'); 
%write cell type on the first line
fprintf(fid, 'GeneSymbol');
numSamp = size(toExp.data,2);
for i = 1:numSamp
    fprintf(fid, '\t');
    fprintf(fid, names{1,i});    
end
fprintf(fid, '\n');    

%write all the data, but always begin with the gene name
numGenes = size(toExp.genes,1);
for g = 1:numGenes
    fprintf(fid, toExp.genes{g,1});
    d = full(toExp.data(g,:));%need to unsparsify for fprintf to work
    for i = 1:numSamp
        fprintf(fid, '\t');
        fprintf(fid, '%f', d(1,i));    
    end
    fprintf(fid, '\n');
end

fclose (fid);

toc




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


