function GenCibersortProfiles(ds, filenameId, nB, nT)

%just assume that the number of cells are more than 100 for both B and T

dsB = ds.cellSubset(ds.cellType == Celltype.BCell);
dsT = ds.cellSubset(ds.cellType == Celltype.TCellCD4Pos | ds.cellType == Celltype.TCellCD8Pos | ds.cellType == Celltype.TCellReg |  ds.cellType == Celltype.TCell);

numGenes = size(dsB.data,1);
numB = size(dsB.data,2);
numT = size(dsT.data,2);

lastB = 1;
lastT = 1;

subData = zeros(numGenes, numB);

for i = 1:nB
    nextB = round(i/nB*numB);
    dsBSub = dsB.cellSubset(lastB:nextB);
    subData(:,i) = mean(dsBSub.data, 2);
end

for i = 1:nT
    nextT = round(i/nT*numT);
    dsTSub = dsT.cellSubset(lastT:nextT);
    subData(:,i + nB) = mean(dsTSub.data, 2);
end

tmpDS = SCDataset();
tmpDS.data = subData;
tmpDS.genes = ds.genes;

WriteDS(tmpDS, strcat(filenameId, '_no_tpm.txt'));
%TPM/CPM normalize:
tmpDS.data = TPM(tmpDS.data);
WriteDS(tmpDS, strcat(filenameId, '.txt'));

    function WriteDS(dsW, filename)
        %Annoyingly, this doesn't work since the variables in the table gets the
        %same name. Do it manually... It will be slow...
        %tmpDS.saveDataTable(filename);
        fid = fopen(filename,'w');
        %first line
        fprintf(fid,'GeneSymbols');
        for k = 1:nB
            fprintf(fid,'\t%s', 'BCell');
        end
        for k = 1:nT
            fprintf(fid,'\t%s', 'TCell');
        end
        fprintf(fid,'\n');

        for row=1:numGenes
            %gene
            fprintf(fid,'%s', dsW.genes{row});

            for col = 1:(nB+nT)
                fprintf(fid,'\t%d', dsW.data(row,col));
            end

            fprintf(fid,'\n');
        end
        fclose(fid);
    end

end
