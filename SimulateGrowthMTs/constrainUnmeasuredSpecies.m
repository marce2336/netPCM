function [statistics, mutantsClasses_wR, listOutliers, growthMTs, processedData] = constrainUnmeasuredSpecies(flagMTsSet, flagSave, flagPruning)
%##########################################################################
%
% With this script the means for the sampled points for all reactions will
% be compared among WT and MTs
%
% flagDataset  -  ('lusk') perform random sampling for Lusk dataset
%                 ('mpi') perform random sampling for MPI dataset
%
% flagSave     -  (0) Default. Don´t save output from processing sampling results
%                 (1) save output from processing sampling results
%
% flagPruning   - (1) Default. Eliminate from lipid profiles the species that are not in model
%                 (0) Keep all species measured in lipidome profiles
%
%##########################################################################

if nargin < 2
    flagSave = 0;
    flagPruning = 1;
    
elseif nargin < 3
    flagPruning = 1;
end

% ###
% ##### 1). Compute means and perform one-sample T-test for WT:
% ###
rootPath = dir(fullfile('..','SamplingResults','Random_samples','tempVariable.mat'));
pathSamples = rootPath.folder;
pathModels = strrep(pathSamples,'Random_samples','Models_sampling');

% load model and samples for WT:
modelSamplingWT = load(fullfile(pathModels,'modelSamplingWT.mat'));
modelSamplingWT = modelSamplingWT.modelSamplingWT;

rootSamplesFile = 'samplesWT_ACHR_';
samplesWT = zeros(numel(modelSamplingWT.rxns),5000);
count = 1;
for i = 1:5
    nameSamplesFile_i = [rootSamplesFile, char(string(i)), '.mat'];
    pathSamples_i = fullfile(pathSamples, nameSamplesFile_i);
    samplesWT_i = load(pathSamples_i);
    samplesWT_i = samplesWT_i.samplesWT_i;
    samplesWT(:, count:count+(size(samplesWT_i,2)-1)) = samplesWT_i;
    count = count+(size(samplesWT_i,2));
end


% compute mean and perform one-sample T-test:
oneSampleTest_WT = zeros(numel(modelSamplingWT.rxns),10);
%labelsCols = {'h', 'p', 'ci_lb', 'ci_ub', 'tstat', 'df', 'sd_lb', 'sd_ub', 'No. Observ', 'meanSamples'};

for i = 1:numel(modelSamplingWT.rxns)
    samplesWT_Rrxni = samplesWT(i,:);
    
    mean_Rxni = mean(samplesWT_Rrxni);
    oneSampleTest_WT(i, end) = mean_Rxni;
    
    oneSampleTest_WT = tTest_oneSample(oneSampleTest_WT,samplesWT_Rrxni,i); 
end





% ###
% ##### 2). Compute means and perform T-tests for MTs:
% ###

switch flagMTsSet
    case {'lusk'}
        % load samples for outlier MTs for Lusk dataset:
        listOutliers = {'AT1G01610';'AT1G02390';'AT1G51260';'AT1G75020';'AT2G19450';'AT2G39290';...
            'AT2G42010';'AT3G18850';'AT4G00400';'AT4G18550';'AT5G03080';'AT5G42870';'AT5G66450'};
     
    case {'mpi'}
        listOutliers = {'AT1G06520';'AT1G73600';'AT1G74320';'AT1G77590';'AT3G02610';'AT3G07690';...
            'AT3G15850';'AT4G00550';'AT4G22340';'AT4G25970';'AT4G27030';'AT4G33030';'AT5G15530';...
            'AT5G19200'};
        
end

        
resultsT_testMTs = '';

for i = 1:numel(listOutliers)
    mutant_i = listOutliers{i};
    
    % load MT model:
    fileModelMT_i = [mutant_i,'_modelSampling_w_relax.mat'];
    modelSamplingMT_i = load(fullfile(pathModels,fileModelMT_i));
    modelSamplingMT_i = modelSamplingMT_i.modelSampling;
    
    % load samples for MT:
    % load sampled fluxes for KO_i:
    newFlag = 'w_relax_';
    samplesMT_i = zeros(numel(modelSamplingMT_i.rxns),5000);
    count = 1;
    for j = 1:5
        nameSamplesFile_j = [mutant_i, '_ACHR_', newFlag, char(string(j)), '.mat'];
        pathSamples_j = fullfile(pathSamples, nameSamplesFile_j);
        samplesMT_j = load(pathSamples_j);
        samplesMT_j = samplesMT_j.sSamples_j;
        samplesMT_i(:, count:count+(size(samplesMT_j,2)-1)) = samplesMT_j;
        count = count+(size(samplesMT_j,2));
    end
    
    resultsT_testMTs.(mutant_i).rxns = modelSamplingMT_i.rxns;
    resultsT_testMTs.(mutant_i).subSystems = modelSamplingMT_i.subSystems;
    
    oneSampleTest_MTi = zeros(numel(modelSamplingWT.rxns),10);
    %labelsCols = {'h', 'p', 'ci_lb', 'ci_ub', 'tstat', 'df', 'sd_lb', 'sd_ub', 'No. Observ', 'meanSamples'};
    twoSampleTest_MTi = zeros(numel(modelSamplingWT.rxns),8);
    %labelsCols = {'h' 'p' 'ci_lb' 'ci_ub' 'tstat' 'df' 'sd_lb' 'sd_ub'};
    
    for j = 1:numel(modelSamplingWT.rxns)
        samplesMTi_Rxnj = samplesMT_i(j,:);
        
        % compute mean of fluxes:
        meanMTi_Rnxj = mean(samplesMTi_Rxnj);
        oneSampleTest_MTi(j, end) = meanMTi_Rnxj;
    
        % carry one-sample T-test for each reaction:
        oneSampleTest_MTi = tTest_oneSample(oneSampleTest_MTi,samplesMTi_Rxnj,j);
        
        % carry out two-sample T-test to find differences among means of WT and MT_i:
        samplesWT_Rxnj = samplesWT(j,:);
        [twoSampleTest_MTi,~] = tTest_twoSample(samplesWT_Rxnj,samplesMTi_Rxnj,twoSampleTest_MTi,j);
    end
    resultsT_testMTs.(mutant_i).oneSampleTtest = oneSampleTest_MTi;
    resultsT_testMTs.(mutant_i).twoSampleTtest = twoSampleTest_MTi;
end



% ###
% ##### 3). Find out for what reaction means are different:
% ###
% Note: since these are the mutants for which no phenotype is expected upon
% the mutation, the results of the two-sample T-test will be used to know 
% in what reactions and subsystems we are getting a differential behavior.

for i = 1:numel(listOutliers)
    mt_i = listOutliers{i};
    results2Ttest_mti = resultsT_testMTs.(mt_i).twoSampleTtest;
    subSystems_mti = resultsT_testMTs.(mt_i).subSystems;
    subSystems_mti = subSystems_mti(1:size(results2Ttest_mti,1));
    
    % get indexes for reactions whose means were different:
    idxDifferentialRxns = results2Ttest_mti(:,1) == 1;
    differentialSubsystems = subSystems_mti(idxDifferentialRxns);
    
    % quantify differential reactions per subsystem:
    uniqueSubsystems = unique(differentialSubsystems);
    hits = zeros(numel(uniqueSubsystems),1);
    for j = 1:numel(uniqueSubsystems)
        subSystem_j = uniqueSubsystems{j};
        idxHits = strcmp(differentialSubsystems, subSystem_j);
        hits(j) = sum(idxHits);
    end
    
    % sort results in descending order:
    [hits,idxs] = sort(hits, 'descend');
    uniqueSubsystems = uniqueSubsystems(idxs);
    resultsT_testMTs.(mt_i).differentialRxns = [uniqueSubsystems,cellstr(string(hits))];
end



% ###
% ##### 4). Insert flux sum for unmeasured lipids:
% ###

% adjust samples for WT by eliminating two samples where ~47% of reactions
% have a zero flux, thus this bias the computation of a threshold
countZerosWT = sum(samplesWT==0);
countZerosWT = find(countZerosWT>1000);
prunedSamplesWT = samplesWT;
prunedSamplesWT(:,countZerosWT) = '';

% get models built for outliers from Lusk dataset and corresponding lipid profiles:
switch flagMTsSet
    case {'lusk'}
        flagMTset='lusk';
        flagDataset=[];
        flagRelaxation=true;
        flagMethod='ACHR';
        flagLipidSet='lusk';
        listMTs = {'AT1G01610';'AT1G02390';'AT1G51260';'AT1G75020';'AT2G19450';'AT2G39290';...
            'AT2G42010';'AT3G18850';'AT4G00400';'AT4G18550';'AT5G03080';'AT5G42870';'AT5G66450'};
        
    case {'mpi'}
        flagMTset='all';
        flagDataset='meanValue_TIC';
        flagRelaxation=true;
        flagMethod='ACHR';
        flagLipidSet='mpi';
        listMTs = {'AT1G06520';'AT1G73600';'AT1G74320';'AT1G77590';'AT3G02610';'AT3G07690';...
            'AT3G15850';'AT4G00550';'AT4G22340';'AT4G25970';'AT4G27030';'AT4G33030';'AT5G15530';...
            'AT5G19200'}; 
end

%[growthMTs_lusk, matchedPhenotypes_lusk, processedData_lusk] = simulateGrowthMT(flagMTset, flagDataset, flagRelaxation, flagMethod, flagLipidSet, listMTs);
[~, matchedPhenotypes, processedData] = simulateGrowthMT(flagMTset, flagDataset, flagRelaxation, flagMethod, flagLipidSet, listMTs, flagPruning);

% load LipidMaps classification:
classificationLM = cellstr(string(readcell(fullfile('InputDataSimulateGrowth', 'LipidMapsClassification.xlsx'))));
codesLM = classificationLM(2:end,4:5);
codesLM = codesLM(:);
codesLM = unique(codesLM);
codesLM = extractBetween(codesLM(2:end), '[', ']');
codesLM = unique(codesLM);

% load reference fluxes for WT:
pFBA_sol = processedData.processedProfiles.phenotypicData.WT.pFBA_sol;

% Pre-allocate space to save predicted growth MTs:
growthMTs = cell(numel(listOutliers)+1,2);
growthMTs(1,:) = {'mutant', 'opt_growth'};

% include data for WT into struct:
mutantData.samplesWT = prunedSamplesWT;

% eliminate from list the outlier MTs for which adding the extra
% constraints, followed by reducing the model, render the model infeasible:
infeasibleMTs = {'AT5G42870','AT2G19450','AT2G39290','AT3G15850'};
idxInfeasible = ismember(listOutliers, infeasibleMTs);
listOutliers(idxInfeasible) = '';
            
for i = 1:numel(listOutliers)
    mt_i = listOutliers{i};
    newName_mti = [mt_i,'outl'];
    growthMTs{i+1,1} = newName_mti;
    
    % for mutant 'AT1G73600' it's necessary to introduce an additional
    % change into the abbreviation for 'MEthP[c]' met before computing sum
    % incoming fluxes:
    switch mt_i
        case {'AT1G73600'}
            modelMT_i = processedData.modelsMT_w.(mt_i);
            idxMet_MTi = strcmp(modelMT_i.mets, 'MEthP[c]');
            modelMT_i.mets(idxMet_MTi) = strrep(modelMT_i.mets(idxMet_MTi), 'MEthP[c]', 'NMEthP[c]');
            processedData.modelsMT_w.(mt_i) = modelMT_i;
            
            idxMet_WTi = strcmp(modelSamplingWT.mets, 'MEthP[c]');
            modelSamplingWT.mets(idxMet_WTi) = strrep(modelSamplingWT.mets(idxMet_WTi), 'MEthP[c]', 'NMEthP[c]');
    end
    
    
    % get sum of incoming fluxes:
    referenceFluxes = getIncomingFluxes(modelSamplingWT, mt_i, codesLM, pFBA_sol);
    
    
    % select lipid pools for which beta values will be obtained (must be determined case-by-case):
    switch mt_i
        case {'AT1G02390','AT3G18850','AT1G51260','AT1G75020','AT4G00400','AT1G01610','AT2G42010',...
                'AT5G23670','AT5G15530','AT4G33030','AT5G16230','AT2G43710','AT3G02610','AT3G15850','AT1G73600','AT1G06520','AT5G19200','AT4G27030','AT3G07690'} 
            % 'AT1G02390' encodes a mitochondrial GPAT isoenzyme that produces 'LPA'
            % Note: in the profiles there are no measurements for 'LPA'
            % species. Here we are going to assume that beta will be equal
            % to the standard deviation of the immediate stream-down
            % product of LPA, that is 'PA'
            % 'AT4G00400' and 'AT1G01610' encode GPAT8 and GPAT 4 catalyzing transference of
            % dicarboxylic FAs into G3P. There's no measurements for these
            % mets. The closest species is DG.
            mutantData.exclude = 'NA';           
            
        case {'AT3G44830','AT3G51520','AT2G19450'}
            % AT3G44830 that catalyzes the synthesis of TG via transference of edited FAs from PC to DG
            % 'AT3G51520' encodes a DGAT enzyme located in the ER
            % 'AT2G19450' encodes DGAT isoenzyme located in ER
            mutantData.exclude = 'LPC';
            mutantData.rxnSuffix = '_r';
            mutantData.poolSuffix = '_TG';
            
        case {'AT5G14180','AT4G18550'} 
            % AT1G54570 encodes PES that degrades TGs under cold stress. This gene is excluded from analysis because during WmPts generation model is infeasible
            % Here we assume that the beta will be equal to the SD of the TGs for which there are measurements available
            % 'AT5G14180' and 'AT4G18550' encodes a lipases producing FFAs. The lipid
            % profiles have no info for FFAs. We will assume sd is equal to the one for the TGs.
            mutantData.exclude = 'LPC';
            
        case {'AT5G66450','AT5G03080','AT1G15080'}
            % 'AT5G66450' and 'AT5G03080' are isoenzymes catalyzing hydrolysis of PA to produce DG in chloroplast.
            % 'AT1G15080' encondes PAP enzyme in chloroplast
            mutantData.exclude = 'NA';
            mutantData.rxnSuffix = '_h';
            mutantData.poolSuffix = '_DG';
                        
        case {'AT5G42870'}
             % 'AT5G42870' is also a PA hydrolase but it's located in the ER.
            mutantData.exclude = 'NA';
            mutantData.rxnSuffix = '_r';
            mutantData.poolSuffix = '_DG';
            
        case {'AT5G57190','AT4G25970'}
            % 'AT5G57190' catalyze conversion of PS into PE.
            %poolID = '] PE ';
            mutantData.exclude = 'NA';
            mutantData.rxnSuffix = '_r';
            mutantData.poolSuffix = '_PE';
            
        case {'AT2G39290'}
            % 'AT2G39290' encodes enzyme converting CDP-DG into PGP in chloroplast
            %poolID = '] PG ';
            mutantData.exclude = 'NA';
            mutantData.rxnSuffix = '_';
            mutantData.poolSuffix = '_PGP';
            
        case {'AT1G54570'}
            mutantData.exclude = 'NA';
            mutantData.rxnSuffix = '_pg';
            mutantData.poolSuffix = '_TG';
            
        case {'AT4G23850'}
            mutantData.exclude = 'NA';
            mutantData.rxnSuffix = '_';
            mutantData.poolSuffix = '_AcylCoA';
            
        case {'AT1G77590'}
            mutantData.exclude = 'NA';
            mutantData.rxnSuffix = '_h';
            mutantData.poolSuffix = '_AcylCoA';
            
        case {'AT4G22340'}
            mutantData.exclude = 'NA';
            mutantData.rxnSuffix = '_r';
            mutantData.poolSuffix = '_CDP_DG';
            
        case {'AT1G74320'}
            mutantData.exclude = 'NA';
            mutantData.rxnSuffix = '_r';
            mutantData.poolSuffix = '_PC';
            
        case {'AT4G00550'}
            mutantData.exclude = 'NA';
            mutantData.rxnSuffix = '_h';
            mutantData.poolSuffix = '_DGDG';
    end
    
    % flagPools  -  ('P') add extra constraints to products
    %               ('R') add extra constraints to reactants
    %               ('B') add extra constraints to products and reactants
    flagPools = 'P';
    
    % add new constraint to selected lipid pool:
    mutantData.newName_mti = newName_mti;
    model_mti = addConstrain2Pool(mt_i,mutantData,processedData,referenceFluxes,flagPools);
    modelsMT.(newName_mti) = model_mti;
    
    % optimize growth of mutant:
    fba_sol = optimizeCbModel_modified(model_mti,'max');
    growthMTs{i+1,2} = cellstr(string(fba_sol.f)); 
    
end
idxEmpty = cellfun(@isempty, growthMTs(:,1));
growthMTs(idxEmpty,:) = '';

% carry-out random sampling:
[feasibleModels,~] = sampleDist4mutants(growthMTs, modelsMT, 'ACHR', true);




% ###
% ##### 5). Process sampling results:
% ###

% get fluxes sampled for MT (both scenarios when available):
fluxesKOs = '';
dirSampledFluxes = dir(fullfile('..','SamplingResults','SampledFluxes','tempVariable.mat'));
path2save = dirSampledFluxes.folder;
files_wRelax = cell(size(growthMTs,1),1);
rootFolder = strrep(path2save, '\SampledFluxes', '');

for i = 2:size(feasibleModels,1)
    ko_i = feasibleModels{i,1};
    locusID = extractBefore(ko_i, 'outl');
    fileKO_i = ['ACHR_wRelax_',ko_i,'.txt'];
    pathKO_i = fullfile(path2save,fileKO_i);
    files_wRelax{i} = ko_i;

    % get list of reactions for WT and retrieve corresponding sampled fluxes:
    idxRxnsWT = contains(modelSamplingWT.grRules, locusID);
    fluxesWT = samplesWT(idxRxnsWT,:);

    % get fluxes sampled for selected KO-rxn with relaxation of bounds:
    sampledFluxesKOi_w = getFluxes4KOrxns(ko_i, fluxesWT, 'w_relax', flagMethod, rootFolder);
    
    % save variable into struct:
    fluxesKOs.(ko_i).wRelaxation = sampledFluxesKOi_w;

    % export results as .txt file:
    writetable(sampledFluxesKOi_w, pathKO_i)

end

path2save_wRelax = fullfile(path2save,'ACHR_wRelax.txt');
idxEmpty = cellfun(@isempty, files_wRelax);
files_wRelax(idxEmpty) = '';
writecell(files_wRelax, path2save_wRelax)




% ###
% ##### STEP 6: Compare the means between WT and MTs:
% ###
dirResults = dir(fullfile('..','SamplingResults','Ttest','tempVariable.mat'));
nameFileStats = ['Results_statistic_tests_Outl_',flagMTsSet,'.xlsx'];
pathStats = fullfile(dirResults.folder, nameFileStats);
listUniqueEnzymes = matchedPhenotypes(strcmp(matchedPhenotypes(:,3), '0'),1);
listIsoenzymes = matchedPhenotypes(strcmp(matchedPhenotypes(:,3), '1'),1);
nameOutlFile = ['Results_classification_mutants_Outl_',flagMTsSet,'.xlsx'];
pathClassifiedMTs = fullfile(dirResults.folder, nameOutlFile);


% carry out T-test for dataset with relaxation, to determine if the
% means are different from zero (0) and if mean_WT is different from mean_MT:
if ~isempty(cellfun(@isempty, files_wRelax))
    listMTs_wRelax = files_wRelax;
    [meanMTs_wR, oneSampleTtest_wR, twoSampleTtest_wR] = perform_Ttest(listMTs_wRelax, 'w_relax',path2save);
    
    statistics.meanMTs_wR = meanMTs_wR;
    statistics.oneSampleTtest_wR = oneSampleTtest_wR;
    statistics.twoSampleTtest_wR = twoSampleTtest_wR;
    writecell(meanMTs_wR, pathStats, 'Sheet', 'meanMTs_wRelaxation')
    writecell(oneSampleTtest_wR, pathStats, 'Sheet', 'oneSampTtest_wRelax')
    writecell(twoSampleTtest_wR, pathStats, 'Sheet', 'twoSampTtest_wRelax')
    
    tsTtestUnique_wR = twoSampleTtest_wR(ismember(twoSampleTtest_wR(:,1), listUniqueEnzymes),:);
    tsTtestUnique_wR = [twoSampleTtest_wR(1,:);tsTtestUnique_wR];
    tsTtestIsoenzymes_wR = twoSampleTtest_wR(ismember(twoSampleTtest_wR(:,1), listIsoenzymes),:);
    tsTtestIsoenzymes_wR = [twoSampleTtest_wR(1,:);tsTtestIsoenzymes_wR];
    
    statistics.tsTtestUnique_wR = tsTtestUnique_wR;
    statistics.tsTtestIsoenzymes_wR = tsTtestIsoenzymes_wR;
    
    % find mutants with full-KO:
    idxFullKO_wR = ~strcmp(oneSampleTtest_wR(:,3), '1');
    idxFullKO_wR(1) = false;
    locusFullKO_wR = oneSampleTtest_wR(idxFullKO_wR,1);
    fullKO_wR = meanMTs_wR(ismember(meanMTs_wR(:,1), locusFullKO_wR),:);
    fullKO_wR = [fullKO_wR, repmat({'full KO'},size(fullKO_wR,1),1)];
    
    % find mutants with full rescue:
    idxFullRescue_wR = strcmp(twoSampleTtest_wR(:,3), '0');
    idxFullRescue_wR(1) = false;
    locusFullRescue_wR = twoSampleTtest_wR(idxFullRescue_wR,1);
    fullRescue_wR = meanMTs_wR(ismember(meanMTs_wR(:,1), locusFullRescue_wR),:);
    fullRescue_wR = [fullRescue_wR, repmat({'full rescue'},size(fullRescue_wR,1),1)];
    
    % find mutants with partial rescue:
    idxPartialRescue_wR = strcmp(twoSampleTtest_wR(:,3), '1');
    idxPartialRescue_wR(1) = false;
    locusPartialRescue_wR = twoSampleTtest_wR(idxPartialRescue_wR,1);
    partialRescue_wR = meanMTs_wR(ismember(meanMTs_wR(:,1), locusPartialRescue_wR),:);
    partialRescue_wR = [partialRescue_wR, repmat({'partial rescue'},size(partialRescue_wR,1),1)];
    
    % create a table with the classification:
    mutantsClasses_wR = [fullKO_wR;fullRescue_wR;partialRescue_wR];
    
    % calculate ratios 'mean-WT/MT':
    ratios_wR = str2double(string(mutantsClasses_wR(:,2)))./str2double(string(mutantsClasses_wR(:,4)));
    mutantsClasses_wR = [mutantsClasses_wR, cellstr(string(ratios_wR))];
    mutantsClasses_wR = [[meanMTs_wR(1,:),{'class' 'mean-WT/MT'}]; mutantsClasses_wR];
    
    mutantsClasses_wR(strcmp(mutantsClasses_wR(:,1),'AT2G29980'),:) = '';
    
    % add phenotypes info:
    mutantClassesLocus = mutantsClasses_wR(:,1);
    mutantClassesLocus = strrep(mutantClassesLocus, 'outl', '');
    [~,idxPhenotype] = ismember(mutantClassesLocus,matchedPhenotypes(:,1));
    idxPhenotype(1) = 1;
    idxPhenotype(idxPhenotype==0) = '';
    mutantsClasses_wR = [mutantsClasses_wR, matchedPhenotypes(idxPhenotype,:)];
    
    
    % save processing results if requested:
    if flagSave == 1
        writecell(mutantsClasses_wR, pathClassifiedMTs, 'Sheet', 'mutantsClasses_wR')
    end
    
end


end