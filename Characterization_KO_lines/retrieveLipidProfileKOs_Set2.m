function processedProfiles = retrieveLipidProfileKOs_Set2(flagPruning)
%##########################################################################
%
% This function filters the raw lipidomics data sets to retrieve the
% profiles for the genes that are included in the merged model
%
% flagFilter - (0) Process lipidomics data for all mutants
%              (1) Default. Process lipidomics data only for mutants of genes in model
%
% flagPruning   - (1) Default. Eliminate from lipid profiles the species that are not in model
%                 (0) Keep all species measured in lipidome profiles
%
%##########################################################################

if nargin < 1
    flagPruning = 1;
end



% ###
% ##### 1). Load lipid profiles
% ###

pathProfiles = fullfile('InputData', 'Lusk_lipids_MTs', 'T_DNA_mutants.xlsx');
lipidProfile = cellstr(string(readcell(pathProfiles, 'Sheet', 'lipid_profiles')));
selectedLipids = readcell(pathProfiles, 'Sheet', 'filtered_lipis');

% filter to select lipid species included in the PLM:
idxSelectedLipid = ismember(lipidProfile(2,:), selectedLipids(2:end,1));
idxSelectedLipid(1:8) = true;
lipidProfile = lipidProfile(:,idxSelectedLipid);
lipidProfile(2,1:8) = lipidProfile(1,1:8);
lipidProfile(1,:) = '';



% ###
% ##### 2). Filter profiles to select genes that are part of Plant LipidModule
% ###
% adjust format for mutants locus ID:
lipidProfile(2:end,5) = strrep(lipidProfile(2:end,5), 't', 'T');
lipidProfile(2:end,5) = strrep(lipidProfile(2:end,5), 'g', 'G');


% get list of genes that are included in the PLM:
dirModel = fullfile('InputData', 'OutputModelUniqueBalanced.mat');
model = load(dirModel);
model = model.model;
modelGenes = model.genes;


% Now find out how many of the KO genes are included in the model:
idxMTs = ismember(lipidProfile(:,5), modelGenes);
idxMTs(1) = true;
lipidsSelectedMTs = lipidProfile(idxMTs, :);
lipidsSelectedMTs = [lipidsSelectedMTs, cell(size(lipidsSelectedMTs,1),1)];
lipidsSelectedMTs{1,size(lipidsSelectedMTs,2)} = 'grRules';

% eliminate mutant 'AT4G11850' for which no measurements were available for any of the replicates:
lipidsSelectedMTs(98:100,:) = '';



% Find out how many genes encode isoenzymes:
for i = 2:size(lipidsSelectedMTs,1)
    gene_i = lipidsSelectedMTs{i,5};
    idxGene_i = contains(model.grRules, gene_i);
    rulesGene_i = unique(model.grRules(idxGene_i));
    lipidsSelectedMTs{i,size(lipidsSelectedMTs,2)} = rulesGene_i{1};
end



% get lipid profiles for Col-0 plants:
idxWT    = contains(lipidProfile(:,5), 'Col-0');
idxWT(1) = true;
lipidsWT = lipidProfile(idxWT,:);



% ###
% ##### 3). Compute averages for biological replicates
% ###
% load information about sum composition:
lipidIDs = selectedLipids;

% filter mutants data:
listDaysMTs  = lipidsSelectedMTs(2:end,2);
locusMTs     = lipidsSelectedMTs(2:end,5);
abundances   = str2double(string(lipidsSelectedMTs(2:end, 9:109)));
uniqueLocus  = unique(locusMTs, 'stable');

% pre-allocate space and perform calculations:
sdRAlipidsMTs       = zeros(numel(uniqueLocus), size(abundances,2)); 
meanIntensityMTs  = sdRAlipidsMTs;
measurementDayMTs = cell(numel(uniqueLocus),1);
bioReplMTs = '';

% calculate average (mean) and standard deviation for MT replicates:
for i = 1:numel(uniqueLocus)
    locusMT_i               = uniqueLocus{i};
    idxLocusMT_i            = strcmp(locusMTs, locusMT_i);
    measurementDayMTs{i}  = unique(listDaysMTs(idxLocusMT_i));

    % calculate mean for replicates
    meanMT_i              = mean(abundances(idxLocusMT_i,:), 1, 'omitnan');
    meanIntensityMTs(i,:) = meanMT_i;

    % calculate SD for relative abundances. When w = 0, the SD is normalized by N-1 (N: number of observations)
    sdLocus_i        = std(abundances(idxLocusMT_i,:), 0, 1, 'omitnan');
    sdRAlipidsMTs(i,:) = sdLocus_i;
    
    % store values for biological replicates:
    bioReplMTs.(locusMT_i).day = unique(listDaysMTs(idxLocusMT_i));
    bioReplMTs.(locusMT_i).MT = abundances(idxLocusMT_i,:)';
end

% calculate averages for Col-0 (wild-type):
listDaysWT   = lipidsWT(2:end,2);
uniqueDaysWT = unique(listDaysWT, 'stable');
abundancesWT = str2double(string(lipidsWT(2:end, 9:109)));
meanWT       = zeros(numel(uniqueDaysWT), size(abundancesWT,2));
bioReplWT = '';

for i = 1:numel(uniqueDaysWT)
    locusWT_i    = uniqueDaysWT{i};
    idxLocusWT_i = strcmp(listDaysWT, locusWT_i);
    meanWT_i     = mean(abundancesWT(idxLocusWT_i,:), 1, 'omitnan');
    meanWT(i,:)  = meanWT_i ;
    
    % store values for biological replicates:
    nameRepWT_i = ['Col0_', locusWT_i];
    bioReplWT.(nameRepWT_i) = abundancesWT(idxLocusWT_i,:)';
end



% ###
% ##### 4). Find out the identity of the lipid species in the model
% ###
pathFx          = fullfile('..','LipidsCoefficientsEstimation','ImplementSLIME','addLipidMapsCodes.m');
pathDir3        = dir(pathFx);
oldDir          = cd(pathDir3.folder);

lipidIDs = [lipidIDs, cell(size(lipidIDs,1),1)];
lipidIDs{1,3} = 'ID model';

for j = 2:size(lipidIDs,1)
    [pos, ~, ~]     = matchToModel(model,lipidIDs{j,2}, true);
    lipidIDs{j,3}   = model.mets(pos);
end

cd(oldDir)
bioReplMTs.lipidIDs = lipidIDs;


% create new lipid labels for data:
listIDs_MTs = lipidsSelectedMTs(1, 9:109);
[~,idxIDs_MTs] = ismember(listIDs_MTs, lipidIDs(:,1));
newListIDs_MTs = lipidIDs(idxIDs_MTs,2:3);




% Eliminate data for lipid species that are not in the model:
switch flagPruning
    case 1
        idxEmpty = cellfun(@isempty, newListIDs_MTs(:,2));
        newListIDs_MTs(idxEmpty,:) = '';
        meanIntensityMTs(:,idxEmpty') = '';
        sdRAlipidsMTs(:,idxEmpty') = '';
        meanWT(:,idxEmpty') = '';
        prunedAbundances = abundances;
        prunedAbundances(:,idxEmpty') = '';

    otherwise
        prunedAbundances = abundances;
        
end



% carry out adjustments for selected species with structural information:
lipids2adjust = {'[GL0201] DG 34:3', 'DAG(16:0_18:3)';
                 '[GL0301] TG 52:5', 'TAG(18:3_34:2)'
                    };


for i = 1:size(lipids2adjust,1)
    species_i = lipids2adjust{i,1};
    splitSpecies_i = strsplit(species_i, ' ');
    idxSpecies_i = contains(newListIDs_MTs(:,1), [splitSpecies_i{1},' ',splitSpecies_i{2},' ',splitSpecies_i{3}]);
    
    switch splitSpecies_i{2}
        case {'DG'} % eliminate species from prokaryotic origin
            modelSpecies_i = newListIDs_MTs{idxSpecies_i, 2};
            modelSpecies_i(contains(modelSpecies_i, 'DGp')) = '';
            newListIDs_MTs{idxSpecies_i, 2} = modelSpecies_i;
            
        case {'TG'} % eliminate species that don't include 18:3-FA
            modelSpecies_i = newListIDs_MTs{idxSpecies_i, 2};
            selectedTGs = true(numel(modelSpecies_i));
            
            for j = 1:numel(modelSpecies_i)
                idxTG_j = strcmp(model.mets, modelSpecies_i{j});
                nameTG_j = model.metNames(idxTG_j);
                
                if ~contains(nameTG_j, '18:3') == 1
                    selectedTGs(j) = false;
                end
            end
            modelSpecies_i(selectedTGs==0) = '';
            newListIDs_MTs{idxSpecies_i, 2} = modelSpecies_i;
    end
end



%%%
%%%%% 3). Calculate the ratios for the relative abundances: RA-MT/RA-WT:
%%%
% Note: the ratios will be calculated using the average of biological replicates
ratiosBioAverages = zeros(numel(measurementDayMTs), size(newListIDs_MTs,1));
ratiosReplicates  = zeros(size(locusMTs,1), size(prunedAbundances,2));
sdFClipidsMTs     = zeros(size(sdRAlipidsMTs,1), size(sdRAlipidsMTs,2));
count = 1;

for i = 1:numel(measurementDayMTs)
    % get measurement for corresponding Col-0 (WT):
    mDay_MTi = measurementDayMTs{i};
    idxWT_i = ismember(uniqueDaysWT, mDay_MTi);
    mIntensityWT_i = mean(meanWT(idxWT_i,:),1);
    
    % calculate ratios for selected MT:
    intensitiesMT_i = meanIntensityMTs(i,:);
    ratiosMT_i = intensitiesMT_i./mIntensityWT_i;
    ratiosBioAverages(i,:) = ratiosMT_i;
    
    % calculate ratios for biological replicates separately:
    locusMT_i = uniqueLocus{i};
    idxLocusMT_i = strcmp(locusMTs, locusMT_i);
    bioReplicates_MTi = prunedAbundances(idxLocusMT_i,:);
    ratiosBioReplicates_MTi = bioReplicates_MTi./mIntensityWT_i;
    ratiosReplicates(count:(count+size(ratiosBioReplicates_MTi,1))-1,:) = ratiosBioReplicates_MTi;
    count = count+size(ratiosBioReplicates_MTi,1);
    
    % calculate SD for FC of relative abundances. When w = 0, the SD is normalized by N-1 (N: number of observations)
    sdFC_MTi            = std(ratiosBioReplicates_MTi, 0, 1, 'omitnan');
    sdFClipidsMTs(i,:)  = sdFC_MTi;
    
    % store values for biological replicates:
    replicatesWT = zeros(numel(listIDs_MTs),24);
    count2 = 1;
    for j = 1:numel(mDay_MTi)
        nameRepWT_i = ['Col0_', mDay_MTi{j}];
        bioReps_j = bioReplWT.(nameRepWT_i);
        replicatesWT(:,count2:count2+size(bioReps_j,2)-1) = bioReps_j;
        count2 = count2+size(bioReps_j,2);
    end
    idxEmpty = sum(replicatesWT==0) == size(replicatesWT,1);
    replicatesWT(:,idxEmpty) = '';
    bioReplMTs.(locusMT_i).WT = replicatesWT;
    
end

% add col and row names to results of biological replicates:
colNames = [{'species'};newListIDs_MTs(:,1)];
FCreplicates = [locusMTs, cellstr(string(ratiosReplicates))];
FCreplicates = [colNames, FCreplicates'];


% compute averages for PI 34:3 because there is data that was acquired in
% positive and negative ion mode:
idxPI34_3 = find(strcmp(newListIDs_MTs(:,1), '[GP0601] PI 34:3'));
newListIDs_MTs(idxPI34_3(2),:) = '';
avMeanPI34_3 = mean(ratiosBioAverages(:,idxPI34_3),2);
ratiosBioAverages(:,idxPI34_3(1)) = avMeanPI34_3;
ratiosBioAverages(:,idxPI34_3(2)) = '';

avSDPI34_3 = mean(sdFClipidsMTs(:,idxPI34_3),2);
sdFClipidsMTs(:,idxPI34_3(1)) = avSDPI34_3;
sdFClipidsMTs(:,idxPI34_3(2)) = '';


% add col and row names to results of averages for biological replicates:
colNames = [{'species'};newListIDs_MTs(:,1)];
FCaverages = [uniqueLocus, cellstr(string(ratiosBioAverages))];
FCaverages = [colNames, FCaverages'];

sdFClipidsMTs = [uniqueLocus, cellstr(string(sdFClipidsMTs))];
sdFClipidsMTs = [colNames, sdFClipidsMTs'];



%%%
%%%%% 4). Save results in struct for further use:
%%%
newListIDs_MTs = [{'species', 'ID model'};newListIDs_MTs];
processedProfiles.meanIntensity  = meanIntensityMTs';
processedProfiles.SD_abundances  = sdRAlipidsMTs';
processedProfiles.lipidIDs       = newListIDs_MTs;
processedProfiles.FCaverages     = FCaverages;
processedProfiles.FCreplicates   = FCreplicates;
processedProfiles.SD_FC          = sdFClipidsMTs;
processedProfiles.KOgenes        = uniqueLocus;
processedProfiles.bioRepsMTs     = bioReplMTs;
processedProfiles.bioRepsWT      = bioReplWT;

end

