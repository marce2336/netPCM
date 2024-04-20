function processedProfiles = processNormLipidProfiles(KOgenes, modelDel, flagSet, flagMeanCol0, flagSave, flagPruning)
%##########################################################################
%
% This function uses as input the normalized lipid profiles and
% assigns the respective codes to the features and eliminates the species
% that are not included in the model
%
% KOgenes       - (cell array) The information for all mutants is provided in this format
%                 ('load') The information must be be created from scratch and loaded!
%
% modelDel      - ('model') If variable is char the model will be created from scratch
%                 (struct)  The model can be provided as struct with COBRA format
%
% flagSet       - ('minValue_median') Dataset processed by filling missing values with 1/5 minimum value 
%                                     for each feature, and normalizing to median.
%
%                 ('meanValue_TIC')   Dataset processed by filling missing values with 1/5 minimum value 
%                                     for each feature, and normalizing to TIC.
%
%                 ('PPCA_median') Dataset processed by filling missing values implementing PPCA 
%                                 and normalizing to median.
%
%                 ('PPCA_TIC')   Dataset processed by filling missing values implementing PPCA 
%                                and normalizing to TIC.
%
% flagMeanCol0  - ('global') Calculates mean of lipids for Col-0
%                           using the intensities for all control plants 
%
%                 ('semiGlobal') calculates the average intensity for the sets of
%                                replicates from Col-0 plants that were measured on the same days
%
%                 ('individual') Default. Calculates the average intensity for each set of
%                                replicates from Col-0 individually, and independently if they were
%                                measured on the same or different days
%
% flagSave      - (true) Default. Save processed lipidomics dataset
%
% flagPruning   - (1) Default. Eliminate from lipid profiles the species that are not in model
%                 (0) Keep all species measured in lipidome profiles
%
%
% PPCA: probabilistic PCA (https://doi.org/10.1093/bioinformatics/btm069)
% TIC: Total Ion Count
%##########################################################################

if nargin < 6
    flagPruning = 1;
end


%%%
%%%%% 1). Obtain normalized lipidomic profiles for Arabidopsis mutants:
%%%

% if information for mutants is not provided, obtain it from scratch:
if  isa(KOgenes, 'char') && isa(modelDel, 'char')
    % Obtain info for mutants:
    KOgenes = findKOmodelGenes(1, 0);

    % obtain model for mutants:
    mutants = load(fullfile('InputData', 'mutants.mat'));
    modelDel = mutants.mutants.modelDel;
end


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

% load running order of samples
runningOrder = cellstr(string(readcell(fullfile('InputData', 'sample_description_RunningOrder.txt'))));

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


% Calculate average (mean) and standard deviation of replicates:
uniqueLabelIDs  = unique(labelIDs(2:end), 'stable');
sdRAlipidsMTs = zeros(size(intensities,1), numel(uniqueLabelIDs)); % pre-allocate space for storing SD
meanIntensity = sdRAlipidsMTs; % pre-allocate space for storing average intensity
setReplicatesMTs = '';
bioReplWT = '';

for i = 1:numel(uniqueLabelIDs)
    label_i             = uniqueLabelIDs{i};
    idxLabel            = strcmp(labelIDs(2:end), label_i);

    % calculate mean for replicates
    meanLabel_i         = mean(intensities(:, idxLabel), 2);
    meanIntensity(:,i)  = meanLabel_i;

    % calculate SD for replicates. When w = 0, the SD is normalized by N-1 (N: number of observations)
    sdLabel_i           = std(intensities(:, idxLabel), 0, 2);
    sdRAlipidsMTs(:,i)  = sdLabel_i;
    
    % store values for biological replicates:
    nameBioRep_i = strrep(label_i, '-','_');
    if contains(label_i, 'KO-')
        setReplicatesMTs.(nameBioRep_i).MT = (intensities(:, idxLabel));
    else
        bioReplWT.(nameBioRep_i).WT = (intensities(:, idxLabel));
    end
    
    
end

% Find out the identity of the lipid species in the model:
pathFx          = fullfile('..','LipidsCoefficientsEstimation','ImplementSLIME','addLipidMapsCodes.m');
pathDir3        = dir(pathFx);
oldDir          = cd(pathDir3.folder);
lipidIDs        = [{'species'}; lipidIDs];
lipidIDs        = addLipidMapsCodes(lipidIDs);
lipidIDs        = [lipidIDs, cell(numel(lipidIDs),1)];
lipidIDs{1,2}   = {'modelID'};

for j = 2:size(lipidIDs,1)
    [pos, ~, ~]     = matchToModel(modelDel,lipidIDs{j}, true);
    lipidIDs{j,2}   = modelDel.mets(pos);
end

cd(oldDir)
bioRepsMTs.lipidIDs = lipidIDs;


% Eliminate species that are not in the model:
switch flagPruning
    case 1
        
        idxEmpty = cellfun(@isempty, lipidIDs(:,2));
        lipidIDs(idxEmpty,:) = '';

        idxEmpty(1) = '';
        meanIntensity(idxEmpty,:) = '';
        intensities(idxEmpty,:) = '';
        sdRAlipidsMTs(idxEmpty,:) = '';
end




%%%
%%%%% 2). Calculate the average of intensities for WT:
%%%

% Here two options can be selected: 
%   (i) ('global') calculates the average intensity considering all samples
%                  and replicates from Col-0.
%
%   (ii) ('semiGlobal') calculates the average intensity for the sets of
%        replicates from Col-0 plants that were measured on the same days
%
%   (iii) ('individual') calculates the average intensity for each set of
%         replicates from Col-0 individually, and independently if they were
%         measured on the same or different days
%
% Note: The averages are later used for the calculation of all ratios:

if ~exist('flagMeanCol0', 'var')
    flagMeanCol0 = 'individual';
end

nameFile = strcat('processedProfiles_', flagSet, '_', flagMeanCol0,'.mat');
namesReplicates = labelIDs(2:end);
switch flagMeanCol0
    case {'global'}
        idxCol0 = contains(namesReplicates, 'Col0');
        meanWT = mean(intensities(:,idxCol0),2);

    case {'semiGlobal'}
        % process information about running order of samples
        runningDays = unique(runningOrder(2:end, 2), 'stable');
        meanWT = zeros(size(intensities,1), numel(runningDays));
        
        for i = 1:numel(runningDays)
            % get day when analysis was carried out:
            runningDay_i = runningDays{i};
        
            % get IDs of Col-0 samples that were run on the selected day:
            idxRunningDay_i = strcmp(runningOrder(:,2), runningDay_i).*contains(runningOrder(:,1), 'Col-0') == 1;
            col0_IDs = runningOrder(idxRunningDay_i,1);
            colNames = mean(str2double(string(strfind(col0_IDs, 'Col-0'))))+4;
            col0_IDs = extractBetween(col0_IDs, 1, colNames);
            col0_IDs = unique(strrep(col0_IDs, '_Col-0', '-Col0'), 'stable');
                
            % calculate mean intensity for selected Col-0 samples:
            idxCol0_i = ismember(namesReplicates, col0_IDs);
            meanCol0_i = mean(intensities(:, idxCol0_i),2);
            meanWT(:,i) = meanCol0_i;
        end

    case {'individual'}
        idxCol0 = contains(uniqueLabelIDs, 'Col0');
        meanWT = meanIntensity(:,idxCol0);
        col0_ids = uniqueLabelIDs(idxCol0);
end



%%%
%%%%% 3). Calculate the ratios for the relative abundances: RA-MT/RA-WT:
%%%
% Note: the ratios will be calculated for the average of replicates
%       ('ratiosAverages') and for each replicate separately
%       ('ratiosReplicates'). 

idxDelReplicates = ~contains(namesReplicates, 'Col0');
labelsReplicates = namesReplicates(idxDelReplicates);

switch flagMeanCol0
    case {'global'} 
        % calculate ratios using global mean RA-WT(Col-0)
        idxAveDel = ~contains(uniqueLabelIDs, 'Col0');
        intensitiesMTs = meanIntensity(:,idxAveDel);
        labelsAverage = uniqueLabelIDs(idxAveDel);
        ratiosAverages = intensitiesMTs./meanWT;
        ratiosReplicates = intensities(:,idxDelReplicates)./meanWT;
        
        % calculate SD for FC of relative abundances. When w = 0, the SD is normalized by N-1 (N: number of observations)
        sdFClipidsMTs     = zeros(size(sdRAlipidsMTs,1), size(sdRAlipidsMTs,2));
        sdFC_MTj            = std(ratiosReplicates, 0, 1, 'omitnan');
        sdFClipidsMTs(i,:)  = sdFC_MTj;

    case {'semiGlobal'} % calculate ratios using mean RA-WT(Col-0) calculated for each running day
        labelsAverage = uniqueLabelIDs(~contains(uniqueLabelIDs, 'Col0'));
        ratiosAverages = zeros(size(meanIntensity,1), numel(labelsAverage));
        ratiosReplicates = zeros(size(intensities,1), sum(idxDelReplicates));
        sdFClipidsMTs     = zeros(size(sdRAlipidsMTs,1), size(sdRAlipidsMTs,2));
        avCount = 1;
        reCount = 1;

        for i = 1:numel(runningDays)
            runningDay_i = runningDays{i};

            % get names of MTs measured each day:
            measuredMTsDay_i = strcmp(runningOrder(:,2), runningDay_i).*(~contains(runningOrder(:,1), 'Col-0')) == 1;
            measuredMTsDay_i = runningOrder(measuredMTsDay_i,1);
            measuredMTsDay_i = unique(extractBefore(measuredMTsDay_i, '_'), 'stable');
            measuredMTsDay_i = strcat('KO-', measuredMTsDay_i);

            % calculate ratios for selected MTs in average intensities:
            idxMTsDay_i = ismember(uniqueLabelIDs, measuredMTsDay_i');
            intensitiesMTs = meanIntensity(:,idxMTsDay_i);
            wtDay_i = meanWT(:,i);
            ratiosDay_i = intensitiesMTs./wtDay_i;
            ratiosAverages(:,avCount:(avCount+size(intensitiesMTs,2))-1) = ratiosDay_i;
            avCount = avCount+size(intensitiesMTs,2);

            % calculate ratios for selected MTs in replicates separately:
            idxMTsDay_i = ismember(namesReplicates, measuredMTsDay_i');
            intensitiesMTs = intensities(:,idxMTsDay_i);
            ratiosDay_i = intensitiesMTs./wtDay_i;
            ratiosReplicates(:,reCount:(reCount+size(intensitiesMTs,2))-1) = ratiosDay_i;
            reCount = reCount+size(intensitiesMTs,2);
            
            % calculate SD for FC of relative abundances. When w = 0, the SD is normalized by N-1 (N: number of observations)
            sdFC_MTj            = std(ratiosDay_i, 0, 1, 'omitnan');
            sdFClipidsMTs(i,:)  = sdFC_MTj;
        end

    case {'individual'} % calculate ratios using RA-WT(Col-0) for each MTs set:
        labelsAverage    = uniqueLabelIDs(idxCol0 == 0);
        ratiosAverages   = zeros(size(meanIntensity,1), numel(labelsAverage));
        sdFClipidsMTs    = ratiosAverages;
        ratiosReplicates = zeros(size(intensities,1), sum(idxDelReplicates));
        lbCount = 1;
        avCount = 1;
        reCount = 1;
        sdCount = 0;

        for i = 1:numel(col0_ids)
            % get subset of mutants for the corresponding control (WT):
            col0_IDi = col0_ids{i};
            findWT_i = strcmp(uniqueLabelIDs, col0_IDi);
            ubCount = find(findWT_i)-1;
            namesMTs = uniqueLabelIDs(lbCount:ubCount);
            
            % calculate ratios for selected MTs in average intensities:
            intensitiesMTs = meanIntensity(:,lbCount:ubCount);
            wtDay_i = meanWT(:,i);
            ratiosDay_i = intensitiesMTs./wtDay_i;
            ratiosAverages(:,avCount:(avCount+size(intensitiesMTs,2))-1) = ratiosDay_i;
            avCount = avCount+size(intensitiesMTs,2);
            lbCount = ubCount+2;

            % calculate ratios for selected MTs in replicates separately:
            idxSelectedMTs = ismember(namesReplicates, namesMTs);
            intensitiesMTs = intensities(:,idxSelectedMTs);
            selectedLocusIDs = namesReplicates(idxSelectedMTs);
            ratiosBioReplicates_MTi = intensitiesMTs./wtDay_i;
            ratiosReplicates(:,reCount:(reCount+size(ratiosBioReplicates_MTi,2))-1) = ratiosBioReplicates_MTi;
            reCount = reCount+size(ratiosBioReplicates_MTi,2);
            
            % calculate SD for FC of relative abundances. When w = 0, the SD is normalized by N-1 (N: number of observations)
            for j = 1:numel(namesMTs)
                nameMT_j                 = namesMTs{j};
                idxMTj                   = strcmp(selectedLocusIDs, nameMT_j);
                subsetRatios_MTj         = ratiosBioReplicates_MTi(:, idxMTj);
                sdFC_MTj                 = std(subsetRatios_MTj, 0, 2, 'omitnan');
                sdCount                  = sdCount + 1;
                sdFClipidsMTs(:,sdCount) = sdFC_MTj;
            end
            
            % store values for biological replicates:
            wti = strrep(col0_IDi, '-','_');
            for k = 1:numel(namesMTs)
                mti = strrep(namesMTs{k}, '-','_');
                setReplicatesMTs.(mti).WT = bioReplWT.(wti).WT;
            end
        
        end      
    
end



%%%
%%%%% 4). Add locus information and adjust MTs with repeated measurements:
%%%

% Assign the corresponding locus ID to the lipid profile of each mutant:
avGeneIDs = addLocusID(labelsAverage, KOgenes);
reGeneIDs = addLocusID(labelsReplicates, KOgenes);

% Double check for samples that were measured twice and calculate mean:
[uniqueIDsRA, locb] = unique(avGeneIDs, 'stable');
FCaverages = zeros(size(ratiosAverages,1), numel(uniqueIDsRA));
sdFCaverages = FCaverages;

for k = 1:numel(uniqueIDsRA)
    gene_i = uniqueIDsRA{k};
    idxReplicates = strcmp(avGeneIDs, gene_i);
    meanAbundance_i = mean(ratiosAverages(:,idxReplicates),2);  
    FCaverages(:,k) = meanAbundance_i;
    
    meanSD_i        = mean(sdFClipidsMTs(:,idxReplicates),2);
    sdFCaverages(:,k) = meanSD_i;
end


% create variable to save ratios calculated using averages:
FCaverages = [uniqueIDsRA; cellstr(string(FCaverages))];
FCaverages = [lipidIDs(:,1), FCaverages];
FCaverages = cellstr(string(FCaverages));

sdFCaverageMTs = [uniqueIDsRA; cellstr(string(sdFCaverages))];
sdFCaverageMTs = [lipidIDs(:,1), sdFCaverageMTs];
sdFCaverageMTs = cellstr(string(sdFCaverageMTs));


% create variable to save ratios calculated for each replicate separately:
FCreplicates = [reGeneIDs; cellstr(string(ratiosReplicates))];
FCreplicates = [lipidIDs, FCreplicates];


% Calculate the CV for mutants that were measured more than once:
getRepeatedLocusMTs = avGeneIDs';
getRepeatedLocusMTs(locb) = '';
cvDuplicatedMTs = cell(size(normIntensities,1)+1, numel(getRepeatedLocusMTs));

for i = 1:numel(getRepeatedLocusMTs)
    locus_i = getRepeatedLocusMTs{i};

    % find id of mutants that were measured twice:
    findReplicates = strcmp(avGeneIDs, locus_i);
    getDuplicatedCodeMTs = labelsAverage(findReplicates);
    findDuplicatedCodeMTs = ismember(labelIDs, getDuplicatedCodeMTs);

    % get the intensities and calculate the CV:
    getDuplicatedintensities = str2double(string(normIntensities(2:end,findDuplicatedCodeMTs)));
    calculateCV = std(getDuplicatedintensities,0,2)./mean(getDuplicatedintensities,2);

    % save results in variable:
    cvDuplicatedMTs{1,i} = locus_i;
    cvDuplicatedMTs{2,i} = join(getDuplicatedCodeMTs, ',');
    cvDuplicatedMTs(3:end,i) = cellstr(string(calculateCV));
end
cvLabels = [normIntensities(1,1); {'duplicatedSamples'}; normIntensities(2:end,1)];
cvDuplicatedMTs = [cvLabels, cvDuplicatedMTs];


% assign locus IDs to uniqueLabelIDs list:
locusUniqueLabelIDs = addLocusID(uniqueLabelIDs, KOgenes);


% assign locus information to biological replicates:
hitsCol0 = ~strcmp(locusUniqueLabelIDs, 'NA');
locusBioReps = [uniqueLabelIDs(hitsCol0); locusUniqueLabelIDs(hitsCol0)];

for i = 1:size(locusBioReps,2)
    fieldsBioRepsMTs = fieldnames(bioRepsMTs);
    koi = strrep(locusBioReps{1,i},'-','_');
    locusi = locusBioReps{2,i};
    idxLocus = strcmp(fieldsBioRepsMTs, locusi);
    
    if sum(idxLocus) > 0
        bioRepsMTs.(locusi).MT = [bioRepsMTs.(locusi).MT, setReplicatesMTs.(koi).MT];
        bioRepsMTs.(locusi).WT = [bioRepsMTs.(locusi).WT, setReplicatesMTs.(koi).WT];
    else
        bioRepsMTs.(locusi).MT = setReplicatesMTs.(koi).MT;
        bioRepsMTs.(locusi).WT = setReplicatesMTs.(koi).WT;
    end
end



%%%
%%%%% 5). Save results in struct for further use:
%%%

processedProfiles.uniqueLabelIDs  = [uniqueLabelIDs;locusUniqueLabelIDs];
processedProfiles.meanIntensity   = meanIntensity;
processedProfiles.SD_abundances   = sdRAlipidsMTs;
processedProfiles.lipidIDs        = lipidIDs;
processedProfiles.FCaverages      = FCaverages;
processedProfiles.FCreplicates    = FCreplicates;
processedProfiles.SD_FC           = sdFCaverageMTs;
processedProfiles.KOgenes         = KOgenes;
processedProfiles.cvDuplicatedMTs = cvDuplicatedMTs;
processedProfiles.bioRepsMTs     = bioRepsMTs;
processedProfiles.bioRepsWT      = bioReplWT;

if ~exist('flagSave', 'var')
    flagSave = true;
end

switch flagSave
    case true
        savePath = fullfile('OutputData', nameFile);
        save(savePath, 'processedProfiles')
end


end
