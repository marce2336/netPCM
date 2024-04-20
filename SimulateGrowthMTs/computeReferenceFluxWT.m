function [phenotypicData, referenceFluxes, processedProfiles] = computeReferenceFluxWT(flagMTset, flagDataset, flagSolver, flagLipidSet, flagPruning)
%##########################################################################
%
% This function calculates the reference flux for the Arabidopsis wild-type
% (Col-0) plants via pFBA. For this, first it is estimated the optimal
% relative growth rate (max v_bio) under the constraints of (i)
% steady-state; (ii) lower and upper flux capacities; (ii) biomass reaction
% under standard conditions (v_(bio,Col0)), (vi) the bound on the ratio
% between the carboxylation and oxygenation reactions catalyzed by RuBisCO
% was set to 2.8876, and (v) the stoichiometric coefficients for chain- and
% backbone-SLIME pseudo reactions were assigned by implementing a custom
% script and using experimemtal information from Col-0 plants under
% standard conditions. Next, FBA was run maximizing growth (biomass rxn was
% set as objective function). This was followed by converting the model
% into irreversible, constraining the lb of the biomass rxn with the
% determined optimal relative growth and computing a reference flux
% distribution (V_Col0) via pFBA. 
% Subsequently, the reactions associated to each lipid species measured for
% the KOs are identified, and the reference fluxes for each of those
% reactions is obtained from the reference flux distribution computed
% above.
%
% flagMTset   -  ('all') load information for the complete set of mutants
%                ('filter') load information only for mutants
%                           corresponding to genes included in the PLM
%                ('lusk') load information for Arabidopsis mutants published by Lusk et al.
%                          (doi:https://doi.org/10.1093/pcp/pcac088)
%
% flagDataset - ('meanValue_TIC')   Dataset processed by filling missing values with 1/5 minimum value 
%                               for each feature, and normalizing to TIC.
%
% flagSolver - ('gurobi')
%              ('ibm_cplex') Default
%
% flagLipidSet - ('mpi') Default. Set of lipid profiles from mutants characterized in-house
%                ('lusk') Set of lipid profiles from Arabidopsis T-DNA lines 
%                         published by Lusk et al. 2022 (Plant Cell Physiol. 63(9): 1193–1204 (2022) 
%                         doi:https://doi.org/10.1093/pcp/pcac088)
%
% flagPruning   - (1) Default. Eliminate from lipid profiles the species that are not included in the Plant Lipid Module
%                 (0) Keep all species measured in lipidome profiles
%
%##########################################################################

if nargin < 5
    flagPruning = 1;
end

% ###
% ##### STEP 1: Retrieve AraCore model expanded with the lipid module:
% ###
% Note: The model was created for Arabidopsis Col-0 under standard conditions
accessionData.Accession         = 1;
accessionData.Day               = 0;
nameDir = fullfile('..','LipidsCoefficientsEstimation','CreateSLIMEmodel','CreateSLIMEModel.m');
pathDir = dir(nameDir);
oldDir  = cd(pathDir.folder);
[biomassRxn, model, ~, ~]    = createSLIMEModel(accessionData,1,1);
phenotypicData.WT.biomassRxn = biomassRxn;



% ###
% ##### STEP 2: Split irreversible reactions into forward and backward reactions:
% ###

% set 'gurobi' as default solver if not specified:
if nargin < 3
    solver = 'cplexlp';
else
    solver = flagSolver;
end
changeCobraSolver(solver, 'LP', 0);


% find optimum biomass and constrain biomass reaction:
s = optimizeCbModel(model);
model.ub(model.c~=0) = s.f;

pathSaveRedModelWT = strrep(oldDir, 'SimulateGrowthMTs', 'SamplingResults');
pathSaveRedModelWT = fullfile(pathSaveRedModelWT,'ReducedModel','reducedModelWT.mat');
dirRedModel = dir(pathSaveRedModelWT);

if isempty(dirRedModel)
    % reduce model:
    [modelRed, hasFlux, maxes, mins] = reduceModel(model);
    
    % obtain irreversible model by running the COBRA function 'convertToIrreversible':
    modelIrrev = convertToIrreversible(modelRed);
    
    % adjust reaction bounds to default values:
    modelIrrev.lb(:) = 0;
    modelIrrev.ub(:) = 1000;

    % find biomass reaction and set as objective function: 
    modelIrrev.c(:) = 0;
    findBiomassRxn = strcmp(modelIrrev.rxns, biomassRxn);
    modelIrrev.c(findBiomassRxn) = 1;

    % save resulst as struct:
    reducedModelWT.model = model;
    reducedModelWT.modelRed = modelRed;
    reducedModelWT.hasFlux = hasFlux;
    reducedModelWT.maxes = maxes;
    reducedModelWT.mins = mins;
    reducedModelWT.modelIrrev = modelIrrev;
    save(pathSaveRedModelWT, 'reducedModelWT')
else
    reducedModelWT = load(pathSaveRedModelWT);
    modelIrrev = reducedModelWT.reducedModelWT.modelIrrev;
end



% ###
% ##### STEP 3: Compute reference flux distribution via pFBA:
% ###   

% calculate Arabidopsis Col-0 flux distribution under control condition via
% pFBA, by minimizing flux through gene associated reactions at optimum
% biomass:
cd(oldDir)
[pFBA_sol, bio_opt] = getFluxDist_pFBA(modelIrrev);
phenotypicData.WT.bio_opt   = bio_opt;
phenotypicData.WT.pFBA_sol  = pFBA_sol;
phenotypicData.model        = model;
phenotypicData.modelIrrev   = modelIrrev;
phenotypicData.bioRxn       = biomassRxn;



% ###
% ##### STEP 4: Get relative abundances of lipids for selected mutants:
% ###
% get list of mutants and corresponding phenotype info:
if nargin < 4
    flagLipidSet = 'mpi';
end

[phenoInfo, listKOs] = getPhenotypeInfo(flagMTset);
phenotypicData.MT.phenoInfo = phenoInfo;
phenotypicData.MT.listKOs = listKOs;
        
nameScript = fullfile('..','Characterization_KO_lines','processNormLipidProfiles.m');
pathScript = dir(nameScript);
oldDir  = cd(pathScript.folder);
        
switch flagLipidSet
    case {'lusk'}
        processedProfiles = retrieveLipidProfileKOs_Set2(flagPruning);
        
    case {'mpi'}
        % get lipid profiles for mutants:
        % Notes: the lipid profiles were (i) pre-processed by filling missing values
        % with 1/5 minimum value for each feature or implementing PPCA; and (ii)
        % normalized to median or total ion count (TIC). The process is described in
        % 'Procedure_normalization_raw_data.txt', located in '..\Paper 5\KO_lines\DataNormalization'
        % The ratios RA-MT/RA-WT was obtained by running the script:
        processedProfiles = processNormLipidProfiles('load', 'model', flagDataset, 'individual', true, flagPruning);

end
cd(oldDir)


% ###
% ##### STEP 5: Find reactions associated to lipid species measured:
% ###

referenceFluxes = '';
for i = 2:size(processedProfiles.lipidIDs,1)
    speciesName = processedProfiles.lipidIDs{i,1};
    fieldName_i = extractAfter(speciesName, '] ');
    fieldName_i = strrep(fieldName_i, ':', '_');
    fieldName_i = strrep(fieldName_i, ' ', '_');
    referenceFluxes.(fieldName_i).speciesID = speciesName;
    speciesInModel = processedProfiles.lipidIDs{i,2};

    for j = 1:numel(speciesInModel)
        species_j = speciesInModel{j};
        
        % find index of each lipid species in model mets list:
        idxSpecies_j = strcmp(modelIrrev.mets,species_j);

        % get stoichiometric coefficients of rxns associated to each lipid species:
        sSpecies_j = full(modelIrrev.S(idxSpecies_j,:));

        % find reactions producing the lipid species_j:
        incomingRxns = sSpecies_j > 0;

        % get reference flux for incoming reactions:
        incomingFluxes = abs(pFBA_sol.v(incomingRxns));

        % save reactions and corresponding fluxes in struct:
        fieldName_j = strrep(species_j, '-', '_');
        fieldName_j = strrep(fieldName_j, '[', '_');
        fieldName_j = strrep(fieldName_j, ']', '');
        referenceFluxes.(fieldName_i).(fieldName_j).speciesID = species_j;
        referenceFluxes.(fieldName_i).(fieldName_j).idxSspecies = find(idxSpecies_j);
        referenceFluxes.(fieldName_i).(fieldName_j).incomingRxns = modelIrrev.rxns(incomingRxns);
        referenceFluxes.(fieldName_i).(fieldName_j).idxIncomingRxns = find(incomingRxns);
        referenceFluxes.(fieldName_i).(fieldName_j).incomingFluxes = incomingFluxes;
    end
end

processedProfiles.phenotypicData = phenotypicData;
        
end

