
dataFolderEF = 'E:/BulkProfiles/BulkProfiles - Copy';

sourceFolder = 'C:/Work/R/RNASeqCTProfileEval';

a = load(strcat(dataFolderEF, '/tmms.mat'));

tmms = Samples();
tmms.data = a.tmmMat;
tmms.sampleIds = a.tmmSampIds;
tmms.genes = a.tmmGenes;

b = load(strcat(dataFolderEF, '/tpms.mat'));
tpms = Samples();
tpms.data = b.tpmMat;
tpms.sampleIds = b.tpmSampIds;
tpms.genes = b.tpmGenes;

c = load(strcat(dataFolderEF, '/qns.mat'));
qns = Samples();
qns.data = c.qnsMat;
qns.sampleIds = c.qnsSampIds;
qns.genes = c.qnsGenes;



%get the design matrix:
dm = xlsread(strcat(sourceFolder, '/DesignMatrix.xlsx'),'DesignMatrix','C3:BX13');





%{
figure
[b,edges] = histcounts(tmms.data(:,1));
xes = zeros(1, size(edges,2)-1);
for i = 1:(size(edges,2)-1)
    xes(1,i) = (edges(1,i) + edges(1,i+1))/2;
end

plot(xes,b)

%}

legs = { 'OC macr.','LIVC T cells', 'PBMC68k T cells', 'LC tumor T cells', 'LC healthy tiss. T cells', 'TCD8 T cells'};

linespec = cell(1,size(dm,2));
%cell type
for i = 1:size(dm,2)
    if dm(9,i) == 2
        linespec{1,i} = 'r';
    else
        linespec{1,i} = 'b';
    end
end

%lab
for i = 1:size(dm,2)
    if dm(2,i) == 1
        linespec{1,i} = 'r';
    elseif dm(2,i) == 2
        linespec{1,i} = 'b';
    elseif dm(2,i) == 3
        linespec{1,i} = 'g';
    elseif dm(2,i) == 5
        linespec{1,i} = 'c';
    elseif dm(2,i) == 6
        linespec{1,i} = 'm';
    end
end


figure
for i = 1:size(tpms.data, 2)
    avg = tpms.data(:,i);
    sel = avg >= 0.5 & avg <= 4000;

    d = log10(avg +0.05);
    a = linspace(-0.3,4,100);

    res = zeros(1, 100);
    for x = 1:size(a,2)
        s = d >= a(x)-0.1 & d <= a(x)+0.1;
        res(1,x) = sum(s);
    end
    res = res ./.2 ./ sum(sel);

    b = 10 .^a;

%    linespec = 'b';
%    if coloring(1,i)
%        linespec = 'r';
%    end
    
    semilogx(b.',res, linespec{1,i});
    hold on;
end

%legend(legs);
xlabel('Gene expression (TPM)')
ylabel('Gene density')
ttl = 'Gene Density';
title(ttl);
axis([0.5 4000 0 0.8]);
set(gca,'FontSize',11);


%TMMs instead
figure
for i = 1:size(tmms.data, 2)
    avg = tmms.data(:,i);
    sel = avg >= 0.5 & avg <= 4000;

    d = log10(avg +0.05);
    a = linspace(-0.3,4,100);

    res = zeros(1, 100);
    for x = 1:size(a,2)
        s = d >= a(x)-0.1 & d <= a(x)+0.1;
        res(1,x) = sum(s);
    end
    res = res ./.2 ./ sum(sel);

    b = 10 .^a;

%    linespec = 'b';
%    if coloring(1,i)
%        linespec = 'r';
%    end
    
    semilogx(b.',res, linespec{1,i});
    hold on;
end

%legend(legs);
xlabel('Gene expression (TMM normalized)')
ylabel('Gene density')
ttl = 'Gene Density';
title(ttl);
axis([0.5 4000 0 0.8]);
set(gca,'FontSize',11);

%quantile normalized data (all should be the same, on top of each other!)
figure
for i = 1:size(qns.data, 2)
    avg = qns.data(:,i);
    sel = avg >= 0.5 & avg <= 4000;

    d = log10(avg +0.05);
    a = linspace(-0.3,4,100);

    res = zeros(1, 100);
    for x = 1:size(a,2)
        s = d >= a(x)-0.1 & d <= a(x)+0.1;
        res(1,x) = sum(s);
    end
    res = res ./.2 ./ sum(sel);

    b = 10 .^a;

    semilogx(b.',res, linespec{1,i});
    hold on;
end

%legend(legs);
xlabel('Gene expression (TPM)')
ylabel('Gene density')
ttl = 'Gene Density';
title(ttl);
axis([0.5 4000 0 0.8]);
set(gca,'FontSize',11);




logdata = LogTrans(tmms.data, 1);

[coeff, scores] = pca(logdata');
figure
gscatter(scores(:,1), scores(:,2), dm(2,:).');
figure
gscatter(scores(:,1), scores(:,2), dm(9,:).');

%technical replicates
figure
techn = dm(11,:).';
techn(isnan(techn)) = 0;
gscatter(scores(:,1), scores(:,2), techn);

%sub cell type
figure
gscatter(scores(:,1), scores(:,2), dm(7,:).');

%tissue
figure
gscatter(scores(:,1), scores(:,2), dm(3,:).');

logdataqn =  LogTrans(qns.data, 1);
[coeffqn, scoresqn] = pca(logdataqn');
figure
gscatter(scoresqn(:,1), scoresqn(:,2), dm(2,:).');
figure
gscatter(scoresqn(:,1), scoresqn(:,2), dm(9,:).');


%some checks of CD8A expression in different cell types
lc = DsLC.get();
bcells = lc.cellSubset(lc.cellType == Celltype.BCell);
macrophages = lc.cellSubset(lc.cellType == Celltype.Macrophage);
tcd8 = lc.cellSubset(lc.cellType == Celltype.TCellCD8Pos);
nk = lc.cellSubset(lc.cellType == Celltype.NKCell);

bcellsMeans = mean(TPM(bcells.data),2);
bcellsMeans(strcmp(bcells.genes, 'CD8A'), :)

macroMeans = mean(TPM(macrophages.data),2);
macroMeans(strcmp(macrophages.genes, 'CD8A'), :)

tcd8Means = mean(TPM(tcd8.data),2);
tcd8Means(strcmp(tcd8.genes, 'CD8A'), :)

nkMeans = mean(TPM(nk.data),2);
nkMeans(strcmp(nk.genes, 'CD8A'), :)




