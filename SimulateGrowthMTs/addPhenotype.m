function matchedPhenotypes = addPhenotype(growthMTs, modelIrrev)
%##########################################################################
%
% This function identifies which of the KO genes encode isoenzymes, and
% make a comparison among the simulated MT growth and the phenotype
% description of mutants collected from different sources
%
% USAGE: matchedPhenotypes = addPhenotype(growthMTs, modelIrrev, phenotypicData)
%
%##########################################################################


% ###
% ##### STEP 1: Identify wich genes code for isoenzymes
% ### 
isIso = repmat({'0'},size(growthMTs,1),1);
isIso{1} = 'is isoenzyme';

for i = 2:size(growthMTs)
    mtLocus_i = growthMTs{i,1};
    getGRrules = (contains(modelIrrev.grRules,mtLocus_i)).*(contains(modelIrrev.grRules, 'or'));

    if sum(getGRrules)>0
        isIso{i} = '1';
    end
end
growthMTs = [growthMTs, isIso];

% eliminate MTs for which lipid profiles were not measured:
idxEmpty = cellfun(@isempty, growthMTs(:,2));
growthMTs(idxEmpty,:) = '';



% ###
% ##### STEP 2: Compare simulated MT growth with phenotype description of mutants
% ### 

% get phenotype information manually curated:
pathFile = fullfile('InputDataSimulateGrowth','Phenotype_selected_KO_genes.xlsx');
phenotypeInfo = cellstr(string(readcell(pathFile, 'Sheet', 'Curated_phenotype_info', 'Range','A1:I78')));
[~, uniqueLocus] = unique(phenotypeInfo(:,1), 'stable');
uniquePhenotype = phenotypeInfo(uniqueLocus,:);

% find for how many MTs there is complete phenotypic information:
[getIdxMTs, loca] = ismember(growthMTs(:,1), uniquePhenotype(:,1));
loca(loca == 0) = '';
matchedPhenotypes = repmat({'NA'},size(growthMTs,1), size(growthMTs,2)+size(uniquePhenotype,2));
matchedPhenotypes(:,1:size(growthMTs,2)) = growthMTs(:,:);
matchedPhenotypes(1,size(growthMTs,2)+1:end) = uniquePhenotype(1,:);
matchedPhenotypes(getIdxMTs,size(growthMTs,2)+1:end) = uniquePhenotype(loca,:);

% identify which of the mutated genes are included in the model:
findModelGenes = ismember(matchedPhenotypes(:,1), modelIrrev.genes);
findModelGenes = cellstr(string(findModelGenes));
findModelGenes{1} = 'ModelGenes';
matchedPhenotypes = [matchedPhenotypes, findModelGenes];
matchedPhenotypes = cellstr(string(matchedPhenotypes));


end

