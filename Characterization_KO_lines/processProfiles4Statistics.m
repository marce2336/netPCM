function [filteredIntensities, filteredGeneIDs] = processProfiles4Statistics(KOgenes, flagSet, growthMTs)
%##########################################################################
%
% This function uses as input the lipid profiles measured in-house for Arabidopsis
% mutant lines. The data is processed by adding the corresponding locus
% IDs.
%
% INPUT:
%    KOgenes   - list T-DNA lines provided as cell array or as char
%
%    flagSet   - ('minValue_median') Dataset processed by filling missing values with 1/5 minimum value 
%                                    for each feature, and normalizing to median.
%   
%                ('meanValue_TIC')   Dataset processed by filling missing values with 1/5 minimum value 
%                                    for each feature, and normalizing to TIC.
%
%                ('PPCA_median') Dataset processed by filling missing values implementing PPCA 
%                                and normalizing to median.
%
%                ('PPCA_TIC')   Dataset processed by filling missing values implementing PPCA 
%                               and normalizing to TIC.
%
%    growthMTs  - cell array with optimal growth values computed for T-DNA lines
%
%
%#########################################################################


%%%
%%%%% 1). Obtain information about Arabidopsis mutants:
%%%
% if information for mutants is not provided, obtain it from scratch:
if  isa(KOgenes, 'char')
    KOgenes = findKOmodelGenes(0, 0);
end


%%%
%%%%% 2). Obtain normalized lipidomic profiles for Arabidopsis mutants:
%%%

% Note: the normalized lipid profiles were obtained as explained in
% 'Procedure_normalization_raw_data.txt', that is located in:
% '..\Characterization_KO_lines\DataNormalization'
switch flagSet
    case {'minValue_median'}
        assignSheetName = 'data_MinValue_normalized_median';  
    case {'meanValue_TIC'}
        assignSheetName = 'data_MinValue_normalized_TIC';
    case {'PPCA_median'}
        assignSheetName = 'data_PPCA_normalized_median';
    case {'PPCA_TIC'}
        assignSheetName = 'data_PPCA_normalized_TIC';
end

pathFile = fullfile('DataNormalization','Data_Normalized','normalized_lipid_profiles.xlsx');
normIntensities = readcell(pathFile, 'Sheet', assignSheetName);

% Add missing labels and harmonize abbreviations usage:
normIntensities{1,1} = 'replicates';
normIntensities      = cellstr(string(normIntensities));
normIntensities(:,1) = strrep(normIntensities(:,1), 'DAG', 'DG');
normIntensities(:,1) = strrep(normIntensities(:,1), 'TAG', 'TG');
normIntensities(:,1) = strrep(normIntensities(:,1), 'Lyso-PC', 'LPC');
normIntensities(:,1) = strrep(normIntensities(:,1), 'Lyso-MGDG', 'MGMG');

% harmonize sample names and get identifier for each group of replicates:
normIntensities(1,:) = strrep(normIntensities(1,:), 'X', '');
normIntensities(1,:) = strrep(normIntensities(1,:), 'Col.0', 'Col-0');

labelIDs = extractBefore(normIntensities(1,:), '_');
findCol0 = contains(normIntensities(1,:), 'Col-0');
labelIDs(findCol0 == 0) = strcat('KO-', labelIDs(findCol0 == 0));
labelIDs(findCol0 == 1) = strcat(labelIDs(findCol0 == 1), '-Col0');
labelIDs{1} = 'group';


% find isomers and add up the intensities:
intensities     = str2double(string(normIntensities(2:end, 2:end)));
findIsomers_1 = contains(normIntensities(2:end,1), '(1)');
findIsomers_2 = contains(normIntensities(2:end,1), '(2)');
addUpIsomers = intensities(findIsomers_1,:) + intensities(findIsomers_2,:);
intensities(findIsomers_1,:) = addUpIsomers;
intensities(findIsomers_2,:) = '';

lipidIDs = normIntensities(2:end,1);
lipidIDs(findIsomers_2) = '';
lipidIDs = strrep(lipidIDs, ' (1)', '');



%%%
%%%%% 3). Add locus information and eliminate mutants that were not simulated:
%%%

% Assign the corresponding locus ID to the lipid profile of each mutant:
replicates = normIntensities(1,:);
labelIDs = [labelIDs;replicates];
geneIDs = addLocusID(labelIDs(:,2:end), KOgenes, 1);

% eliminate mutants that were not included in the simulations:
if nargin == 3
    idxSimulatedMTs = ismember(geneIDs(2,:), growthMTs(:,1));
    filteredIntensities = intensities(:,idxSimulatedMTs);
    filteredGeneIDs = geneIDs(:,idxSimulatedMTs);
else
    filteredGeneIDs = geneIDs(2,:);
    idxSimulatedMTs = true(numel(filteredGeneIDs));
    filteredIntensities = intensities;
    filteredGeneIDs = geneIDs;
end


% select control samples (Col-0) corresponding to the selected species:
idxCol0 = find(contains(geneIDs(2,:), 'Col0'));
setsCol0 = zeros(1,numel(idxCol0));
for i = 1:(numel(idxCol0)-1)
    setsCol0(i) = idxCol0(i+1) - idxCol0(i);
end

% find out which samples correspond to standard (Col-0) samples:
setsCol0 = setsCol0>1;
idxSetsCol0 = zeros(sum(idxCol0),2);
count = 0;
for i = 1:numel(setsCol0)
    idx_i = setsCol0(i);
    if i == 1
        idxSetsCol0(1,1) = idxCol0(i);
        
    else
        if idx_i == true
            count = count+1;
            idxSetsCol0(count,2) = idxCol0(i);
            idxSetsCol0(count+1,1) = idxCol0(i+1);
        end
    end
end
idxSetsCol0(count+1,2) = idxCol0(i);
idxSetsCol0(sum(idxSetsCol0,2)==0,:) = '';


% identify which standards to leave according to mutants simulated:
col0ToKeep = cell(size(idxSetsCol0,1),1);
count = 1;
for i = 1:size(idxSetsCol0,1)
    if i == 1
        simulatedMTs = sum(idxSimulatedMTs(1:(idxSetsCol0(i,1)-1)));
        if simulatedMTs > 1
            col0ToKeep{i} = geneIDs(2,idxSetsCol0(i,1));
        end
    else
        simulatedMTs = sum(idxSimulatedMTs((idxSetsCol0(i-1,2)+1):(idxSetsCol0(i,1)-1)));
        if simulatedMTs > 1
            count = count+1;
            col0ToKeep{count} = geneIDs(2,idxSetsCol0(i,1));
        end
    end
end
idxEmpty = ~cellfun(@isempty, col0ToKeep);
col0ToKeep = cellstr(string(col0ToKeep(idxEmpty)));

% calculate average for standards:
meanStandards = zeros(size(filteredIntensities,1), numel(col0ToKeep));
for i = 1:numel(col0ToKeep)
    col0_i = col0ToKeep{i};
    idxCol0_i = strcmp(geneIDs(2,:), col0_i);
    meanCol0_i = mean(intensities(:,idxCol0_i),2);
    meanStandards(:,i) = meanCol0_i;
end
col0Labels = extractBefore(col0ToKeep, '-Col0');
col0Labels = strrep(col0Labels, 'T', 'Col0_');
col0Labels = [col0Labels'; repmat({'Col-0'},1,numel(col0Labels))];

% group mutants for which more than one group set was measured:
[uniqueLocusIDs, ~] = unique(filteredGeneIDs(2,:), 'stable');
for i = 1:numel(uniqueLocusIDs)
    uLocus_i = uniqueLocusIDs{i};
    idxLocus_i = strcmp(filteredGeneIDs(2,:), uLocus_i);
    suffixLocus_i = cellstr(string(1:sum(idxLocus_i)));
    suffixLocus_i = strcat(uLocus_i,'_',suffixLocus_i);
    filteredGeneIDs(1,idxLocus_i) = suffixLocus_i;
end

% create one table with results and save as .csv file:
filteredIntensities = [lipidIDs, cellstr(string(filteredIntensities)), cellstr(string(meanStandards))];
filteredGeneIDs = [filteredGeneIDs, col0Labels];


end

