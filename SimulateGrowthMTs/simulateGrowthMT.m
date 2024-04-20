function [growthMTs, matchedPhenotypes, processedData] = simulateGrowthMT(flagMTset, flagDataset, flagRelaxation, flagMethod, flagLipidSet, listMTs, flagPruning)
%##########################################################################
% 
% USAGE:
%   [growthMTs, matchedPhenotypes, processedData, fluxesKOs, checkGrowth, statistics] = simulateGrowthMT(flagMTset, flagDataset, flagRelaxation, flagMethod, flagLipidSet, listMTs, flagPruning)
%
% flagMTset   -  ('all') load information for the complete set of mutants
%                ('filter') load information only for mutants
%                           corresponding to genes included in the Plant Lipid Module
%                ('lusk') load information for Arabidopsis mutants published by Lusk et al.
%                          (doi:https://doi.org/10.1093/pcp/pcac088)
%
% flagDataset - ('meanValue_TIC')   Dataset processed by filling missing values with 1/5 minimum value 
%                           for each feature, and normalizing to total ion count (TIC)              
%
% flagRelaxation - (true) Default. Use standard deviation (sd) to relax lower and upper bounds
%                  (false) Do not relax upper and lower bounds
%
% flagMethod - ('ACHR') Default. Use Artificial Centered Hit-and-Run (ACHR) included in the CobraToolbox  
%
% flagLipidSet - ('mpi') Default. Set of lipid profiles from mutants characterized in-house
%                ('lusk') Set of lipid profiles from Arabidopsis T-DNA lines 
%                         published by Lusk et al. 2022 (Plant Cell Physiol. 63(9): 1193–1204 (2022) 
%                         doi:https://doi.org/10.1093/pcp/pcac088)
%
% flagPruning   - (1) Default. Eliminate from lipid profiles the species that are not included in the Plant Lipid Module
%                 (0) Keep all species listed in lipidome profiles
%
%##########################################################################

% ###
% ##### STEP 1: Compute reference flux distribution for WT:
% ###
% NOTE: the default method is 'meanValue_TIC'.
if nargin == 0
    flagMTset = 'all';
    flagDataset = 'meanValue_TIC';
    flagRelaxation = true;
    flagMethod = 'ACHR';
    flagLipidSet = 'mpi';
    flagPruning = 1;
    
elseif nargin < 2
    flagDataset = 'meanValue_TIC';
    flagRelaxation = true;
    flagMethod = 'ACHR';
    flagLipidSet = 'mpi';
    flagPruning = 1;

elseif nargin < 3
    flagRelaxation = true;
    flagMethod = 'ACHR';
    flagLipidSet = 'mpi';
    flagPruning = 1;

elseif nargin < 4
    flagMethod = 'ACHR';
    flagLipidSet = 'mpi';
    flagPruning = 1;
    
elseif nargin < 5
    flagLipidSet = 'mpi';
    flagPruning = 1;
    
elseif nargin < 7
    flagPruning = 1;
end

[phenotypicData, referenceFluxes, processedProfiles] = computeReferenceFluxWT(flagMTset, flagDataset, 'ibm_cplex', flagLipidSet, flagPruning);

modelIrrev = phenotypicData.modelIrrev;

% create struct to store data necessary to compute sum flux for mutants:
processedData.processedProfiles = processedProfiles;
processedData.referenceFluxes = referenceFluxes;
processedData.modelIrrev = modelIrrev;
processedData.flagRelaxation = flagRelaxation;





% ###
% ##### STEP 2: Adjust list of mutants:
% ###

% eliminate duplicated entries in mutants list:
listKOs = phenotypicData.MT.listKOs;
[~, locb] = unique(listKOs(:,2), 'stable');
uniqueListKOs = listKOs(locb,:);

% eliminate genes for which lipid profiles were not measured:
notMeasuredMTs = ~ismember(uniqueListKOs(:,2), processedProfiles.FCaverages(1,:)');
notMeasuredMTs(1) = false;
uniqueListKOs(notMeasuredMTs,:) = '';


% eliminate mutants that correspond to:
%   (i) genes that are not included in the extended model
modelGenes = modelIrrev.genes;
commonGenes = ismember(uniqueListKOs(:,2), modelGenes);
commonGenes(1) = true;
uniqueListKOs = uniqueListKOs(commonGenes,:); 

%   (ii) genes that are not expressed in rossette or
%   (iii) genes whose mutation is embryo lethal in homozygous mutants
%               'locusID'      'expression profile'             'phenotype'
genes2remove = {'AT2G46210'     'not expressed in rossette'     'Growth affected mostly in roots';
                'AT3G13930'     'expressed in leaves'           'Embryo lethal in homozygous mutants';
                'AT3G16950'     'expressed in leaves'           'Embryo lethal in homozygous mutants';
                'AT4G09760'     'expressed in leaves'           'Growth affected mostly in roots';
                'AT4G29460'     'Not expressed in rossette'     'Not expressed in rossette';
                'AT4G30580'     'expressed in leaves'           'Embryo lethal in homozygous mutants'};

idxGenes2remove = ismember(uniqueListKOs(:,2), genes2remove(:,1));
uniqueListKOs(idxGenes2remove,:) = '';





% ###
% ##### STEP 3: Calculate flux sum-MTs, constrain models and optimize growth:
% ###

% Note 1: The case scenario selected for the analysis is (1)
% Note 2: The analysis will be carried out with relaxation of bounds


% Optimize growth relaxing the lower and upper bounds using the sd:
processedData.flagRelaxation = true;
[growthMTs, ~, modelsMT_w] = optimizeGrowthMT(uniqueListKOs, processedData, 1, 'ibm_cplex');
processedData.modelsMT_w = modelsMT_w;
growthMTs = cellstr(string(growthMTs));

% filter list of mutants to simulate if a list is provided:
if nargin > 5 && ~isempty(listMTs)
    idxMTs = ismember(growthMTs(:,1), listMTs');
    idxMTs(1) = true;
    growthMTs = growthMTs(idxMTs,:);
end





% ###
% ##### STEP 4: Compare simulated MT growth with phenotype description of mutants
% ### 
matchedPhenotypes = addPhenotype(growthMTs, modelIrrev);





% ###
% ##### STEP 5: Sampling the solution space for WT:
% ### 

% get flux intervals previously computed for WT
pathFileFVA_WT = fullfile('..', 'SamplingResults', 'FVA_results', 'fva_WT_na.mat');
dirFileFVA_WT = dir(pathFileFVA_WT);
modelWT = modelIrrev;

if isempty(dirFileFVA_WT)% constrain biomass reaction with optimal growth and perform FVA
    modelWT.ub(modelWT.c~=0) = phenotypicData.WT.bio_opt;
    [minFluxWT, maxFluxWT] = fluxVariability(modelWT,100,'max');
    fva_WT.minFlux = minFluxWT;
    fva_WT.maxFlux = maxFluxWT;
    save(pathFileFVA_WT, 'fva_WT')
else % load FVA results already computed
    fva_WT = load(pathFileFVA_WT);
    minFluxWT = fva_WT.fva_WT.minFlux;
    maxFluxWT = fva_WT.fva_WT.maxFlux;
end

% select method for sampling the solution space:
modelWT.ub(modelWT.c~=0) = 1000;

switch flagMethod 
        case {'ACHR'}
            nameExistingFile = 'samplesWT_ACHR_1.mat';
            dirSamples = dir(fullfile('..', 'SamplingResults','Random_samples','tempVariable.mat'));
            pathWTsamples = fullfile(dirSamples.folder, nameExistingFile);
            fileName = extractBefore(nameExistingFile, '_1.mat');


            % run 'sampleCbModel' if the analysis was not carried out:
            if isempty(dir(pathWTsamples))
                options.nStepsPerPoint = 200;
                options.nPointsReturned = 5000;
                options.nWarmupPoints = 15000;
                options.nFiles = 10;
                options.nPointsPerFile = 1000;
                options.nFilesSkipped = 2;
                options.minFlux = minFluxWT;
                options.maxFlux = maxFluxWT;
                options.bioOpt = phenotypicData.WT.bio_opt;
                options.mutant_i = 'WT';
                options.flagRelaxation = 'na';
                options.flagModelReduction = false;
                options.routeFolder = fullfile('OutputData','ACHR','FilesSamplingProcedure','WT');
                [modelSamplingWT, samplesWT, ~] = sampleCbModel_modified(modelWT, fileName, 'ACHR', options);
                
                %split samples to reduce file size and save samples:
                count = 1;
                for i = 1:5
                    samples_i = samplesWT(:,count:(1000*i));
                    nameSamples_i = ['samplesWT_ACHR_', char(string(i)), '.mat'];
                    pathSamples_i = fullfile(dirSamples.folder, nameSamples_i);
                    save(pathSamples_i, 'samples_i')
                    count = count+(1000);
                end
            
                % save the model sampling for later use!
                pathModelSampling = strrep(dirSamples.folder, 'Random_samples', 'Models_sampling');
                pathModelSampling = fullfile(pathModelSampling,'modelSamplingWT.mat');
                save(pathModelSampling, 'modelSamplingWT')
                
            end
end





% ###
% ##### STEP 6: Sampling the solution space for MTs:
% ###

switch flagMethod
    case {'ACHR'}
        [~,~] = sampleDist4mutants(growthMTs, modelsMT_w, flagMethod);
       
end


end

