function statsLipidsOutliers = compareLipidProfiles(outliers_w_lipids,processedData_mpi,processedData_lusk,pathSave)
%##########################################################################
%
% Performs a comparison of the abundances for the lipid species measured in
% the MTs against the WT-plants by mean of a two-sample T-test
%
%##########################################################################

cDir = pwd;
dirTtest = strrep(cDir, 'Figures_source_data\Figure_2_&_Supplemental_Figure_1', 'SimulateGrowthMTs');
oldDir = cd(dirTtest);

statsLipidsOutliers = '';
for i = 1:size(outliers_w_lipids,1)
    mt_i = outliers_w_lipids{i,1};
    dataset_i = outliers_w_lipids{i,2};
    species_i = outliers_w_lipids{i,3};
    
    % get lipid profiles measured for selected mutant:
    switch dataset_i
        case {'mpi'}
            lipidIDs = processedData_mpi.processedProfiles.bioRepsMTs.lipidIDs;
            lipidIDs(contains(lipidIDs(:,1),'PG 34:4')) = strrep(lipidIDs(contains(lipidIDs(:,1),'PG 34:4')), 'PG 34:4 ', 'PG 34:4');
            lipidProfile_mti = processedData_mpi.processedProfiles.bioRepsMTs.(mt_i).MT;
            lipidProfile_wti = processedData_mpi.processedProfiles.bioRepsMTs.(mt_i).WT;
            sumComposition = lipidIDs(2:end,1);
            meanFCs = processedData_mpi.processedProfiles.FCaverages;
            meanFCs(2:end,1) = extractAfter(meanFCs(2:end,1), '] ');
            
            % remove uninteresting species to improve heatmap visualization:
            species2rmv = {'[GL0501] DGDG 32:0','[GL0501] DGDG 32:2','[GL0501] DGDG 32:3',...
                '[GL0501] MGDG 32:2','[GL0501] MGDG 32:3','[GL0501] MGDG 32:4','[GL0501] MGDG 32:5',...
                '[GL0501] MGDG 32:6'};

            idx2rmv = ismember(lipidIDs(:,1), species2rmv);
            lipidIDs(idx2rmv,:) = '';
            meanFCs(idx2rmv,:) = '';
            idx2rmv(1) = '';
            lipidProfile_mti(idx2rmv,:) = '';
            lipidProfile_wti(idx2rmv,:) = '';
            sumComposition(idx2rmv,:) = '';
            
            
        case {'lusk'}
            lipidIDs = processedData_lusk.processedProfiles.bioRepsMTs.lipidIDs;
            lipidProfile_mti = processedData_lusk.processedProfiles.bioRepsMTs.(mt_i).MT;
            lipidProfile_wti = processedData_lusk.processedProfiles.bioRepsMTs.(mt_i).WT;
            sumComposition = lipidIDs(2:end,2);
            meanFCs = processedData_lusk.processedProfiles.FCaverages;
            meanFCs(2:end,1) = extractAfter(meanFCs(2:end,1), '] ');
            
            % carry out manual adjustments to eliminate duplicated PE 42:3 by computing their average:
            idxPE = strcmp(meanFCs(:,1), 'PE 42:3');
            avPE = mean(str2double(string(meanFCs(idxPE,2:end))));
            idxPE = find(idxPE);
            meanFCs(idxPE(1),2:end) = cellstr(string(avPE));
            meanFCs(idxPE(2),:) = '';

    end
    meanFCs(contains(meanFCs(:,1),'PG 34:4')) = strrep(meanFCs(contains(meanFCs(:,1),'PG 34:4')), 'PG 34:4 ', 'PG 34:4');
    
    % compare abundances of MT vs. WT for selected species:
    labelsCols = {'species' 'h' 'p' 'ci_lb' 'ci_ub' 'tstat' 'df' 'sd_lb' 'sd_ub'};
    switch species_i
        case {'all'}
            twoSampleTest_mti = zeros(size(lipidProfile_mti,1),8);
            labels_mti = cell(size(lipidProfile_mti,1),2);
            prunedLipidProfile_mti = lipidProfile_mti;
            prunedLipidProfile_wti = lipidProfile_wti;
            rowNames_mti = sumComposition;
            
        otherwise
            % get list of lipid classes:
            listSpecies_mti = split(sumComposition, ' ');
            idxSpecies_i = ismember(listSpecies_mti(:,2), split(species_i, ','));
            prunedLipidProfile_mti = lipidProfile_mti(idxSpecies_i,:);
            prunedLipidProfile_wti = lipidProfile_wti(idxSpecies_i,:);
            rowNames_mti = sumComposition(idxSpecies_i);
            
            twoSampleTest_mti = zeros(size(prunedLipidProfile_mti,1),8);
            labels_mti = cell(size(prunedLipidProfile_mti,1),2);
    end
    
    
    % perform two-sample T-test:
    for j = 1:size(prunedLipidProfile_mti,1)
        % perform two-sample T-test:
        samplesMTi = prunedLipidProfile_mti(j,:);
        samplesWTi = prunedLipidProfile_wti(j,:);
        [twoSampleTest_mti,~] = tTest_twoSample(samplesWTi,samplesMTi,twoSampleTest_mti, j);

        % get lipid ID and FC computed for RA-MT/RA-WT:
        species_j = extractAfter(rowNames_mti{j}, '] ');
        idxSpecies_j = strcmp(meanFCs(:,1), species_j);
        idxLocus_i = strcmp(meanFCs(1,:), mt_i);
        fc_mti = meanFCs(idxSpecies_j, idxLocus_i);
        labels_mti(j,1) = fc_mti;
        
        if str2double(string(fc_mti)) > 1%fc_mti > 1
            labels_mti{j,2} = 'increase';
        else
            labels_mti{j,2} = 'decrease';
        end
    end
    twoSampleTest_mti = [rowNames_mti,cellstr(string(twoSampleTest_mti))];
    twoSampleTest_mti = [labelsCols;twoSampleTest_mti];
    labels_mti = [{'FC(MT/WT)' 'direction'};labels_mti];
    twoSampleTest_mti = [twoSampleTest_mti,labels_mti];
    
    % save results into struct:
    statsLipidsOutliers.(mt_i).species = species_i;
    statsLipidsOutliers.(mt_i).lipidIDs = lipidIDs;
    statsLipidsOutliers.(mt_i).lipidProfile_mt = lipidProfile_mti;
    statsLipidsOutliers.(mt_i).lipidProfile_wt = lipidProfile_wti;
    statsLipidsOutliers.(mt_i).twoSampleTest = twoSampleTest_mti;
    
    % export results as .xlsx file:
    writecell(twoSampleTest_mti, pathSave, 'Sheet', mt_i)
    
end

cd(oldDir)

end

