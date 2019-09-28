
%% 1. Create single-cell profiles for CIBERSORTx. Use all in the cancer part of the lc dataset
[lct, ~] = DsLC.get();

unique(CelltypeId2CelltypeName(lct.cellType))
figure
histogram(lct.cellType)

%remove rare cell types for which there will not be eonough cells to pool
lct = lct.cellSubset((lct.cellType ~= Celltype.Dendritic) & (lct.cellType ~= Celltype.Epithelial) & (lct.cellType ~= Celltype.Langerhans) );

%The patient cancer types:
%1 Female 70y Yes (2) pT2bN0M0 IIA Squamous Left upper Active
%2 Male 86y Yes (2) pT2bN0M0 IB Squamous Right upper Former
%3 Male 68y No pT4N2M0 IIIB Adenomatous Right upper Former
%4 Female 64y Yes (2) pT2aN1M0 IIB Adenomatous Left upper Former
%5 Male 60y Yes (2) pT1cN0M0 IA3 Large cell Left upper Former
%6 Male 65y Yes (2) pT4N1M0 IIIA Adenomatous Left upper Former
%7 Male 60y No pT2aN0M0 IB Squamous Left upper Former
%8 Female 55y No pT3N0M0 IIB Pleiomorphic Right upper Active
%so, pat 1 and 2 are LUSC, pat 3 and 4 LUAD

%It turns out we have only a small number of malignant cells from LUSC (400), so only 
%work with LUAD malignant cells and remove all the others
selLuad = (lct.cellType ~= Celltype.Malignant) | strcmp(lct.sampleIds,'3') | strcmp(lct.sampleIds,'4');
lct = lct.cellSubset(selLuad);


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



%selLusc = (lct.cellType == Celltype.Malignant) & (strcmp(lct.sampleIds,'1') | strcmp(lct.sampleIds,'2'))
%sum(selLusc)
%sum(selLuad)
%lct.cellType(1,selLusc) = 300;
%lct.cellType(1,selLuad) = 301;
%names(1, lct.cellType == 300) = {'MalLUSC'};
%names(1, lct.cellType == 301) = {'MalLUAD'};


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

%make sure there are no spaces in the names
names = CelltypeId2CelltypeName(toExp.cellType);
names(1, toExp.cellType == Celltype.BCell) = {'BCell'};
names(1, toExp.cellType == Celltype.TCell) = {'TCell'};
names(1, toExp.cellType == Celltype.TCellCD4Pos) = {'TCD4'};
names(1, toExp.cellType == Celltype.TCellCD8Pos) = {'TCD8'};
names(1, toExp.cellType == Celltype.TCellReg) = {'TReg'};
names(1, toExp.cellType == Celltype.NKCell) = {'NKCell'};

unique(names)



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


