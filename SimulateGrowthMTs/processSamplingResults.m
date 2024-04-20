function [statistics, processedData, growthMTs] = processSamplingResults(flagDataset, flagPruning)
%##########################################################################
%
% The script randomly samples the feasible solution space for models
% created for Arabidopsis T-DNA lines. The models are constrained with the 
% use of the relative abundance of measured lipids by using them to
% constrain the ratio of the respective flux sums between the wild type and
% mutant.
%
% flagDataset   -  ('all') Default. load information for the complete set of mutants
%                  ('mpi') load information only for mutants measured at MPI
%                  ('lusk') load information only for mutants published by Lusk et al.
%                          (doi:https://doi.org/10.1093/pcp/pcac088)
%
% flagPruning   - (1) Default. Eliminate from lipid profiles the species that are not in model
%                 (0) Keep all species measured in lipidome profiles
%
%##########################################################################

% ###
% ##### STEP 1: Obtain random samples for Arabidopsis wild type Col-0 and T-DNA lines
% ###

listMTs_mpi = {'AT1G06520','AT1G08510','AT1G63050','AT1G73600','AT1G74320','AT1G74960',...
    'AT1G77590','AT2G32260','AT2G43710','AT3G02610','AT3G07690','AT3G15850','AT3G18000',...
    'AT3G48780','AT4G00550','AT4G22340','AT4G23850','AT4G25970','AT4G27030','AT4G33030',...
    'AT5G15530','AT5G16230','AT5G19200','AT5G23670','AT5G46290'};


listMTs_lusk = {'AT1G01610','AT1G02390','AT1G15080','AT1G32200','AT1G51260','AT1G52760',...
    'AT1G54570','AT1G64670','AT1G75020','AT1G80950','AT2G19450','AT2G37940','AT2G39290',...
    'AT2G42010','AT3G14360','AT3G15730','AT3G18850','AT3G26840','AT3G44830','AT3G51520',...
    'AT3G56940','AT3G57140','AT4G00400','AT4G18550','AT4G26910','AT5G03080','AT5G04040',...
    'AT5G14180','AT5G25370','AT5G37300','AT5G42870','AT5G57190','AT5G66450','AT1G12640',...
    'AT3G05630','AT1G67560','AT1G72520','AT2G38110','AT3G03540','AT3G16785','AT3G57650',...
    'AT3G58490','AT3G62860','AT4G00240','AT4G01950','AT5G14470','AT5G16120'};

    
if nargin < 1
    flagDataset = 'all';
    flagPruning = 1;
    
elseif nargin < 2
    flagPruning = 1;
end

        
switch flagDataset
    case {'all'}
        % simulate growth for MPI dataset:
        flagMTset='all';
        flagDataset='meanValue_TIC';
        flagRelaxation=true;
        flagMethod='ACHR';
        flagLipidSet='mpi';
        [growthMTs_mpi, matchedPhenotypes_mpi, processedData_mpi] = simulateGrowthMT(flagMTset, flagDataset, flagRelaxation, flagMethod, flagLipidSet,listMTs_mpi,flagPruning);

        % simulate growth for Lusk dataset:
        flagMTset='lusk';
        flagDataset=[];
        flagRelaxation=true;
        flagMethod='ACHR';
        flagLipidSet='lusk';
        [growthMTs_lusk, matchedPhenotypes_lusk, processedData_lusk] = simulateGrowthMT(flagMTset, flagDataset, flagRelaxation, flagMethod, flagLipidSet,listMTs_lusk,flagPruning);

        % merge results:
        growthMTs = [growthMTs_mpi; growthMTs_lusk(2:end,:)];
        [~,idxUnique] = unique(growthMTs(:,1), 'stable');
        growthMTs = growthMTs(idxUnique,:);
        
        matchedPhenotypes = [matchedPhenotypes_mpi;matchedPhenotypes_lusk(2:end,:)];
        [~,idxUnique] = unique(matchedPhenotypes(:,1), 'stable');
        matchedPhenotypes = matchedPhenotypes(idxUnique,:);
        
        processedData.setMPI = processedData_mpi;
        processedData.setLusk = processedData_lusk;
        
        
    case {'mpi'}
        flagMTset='all';
        flagDataset='meanValue_TIC';
        flagRelaxation=true;
        flagMethod='ACHR';
        flagLipidSet='mpi';
        [growthMTs_mpi, matchedPhenotypes_mpi, processedData_mpi] = simulateGrowthMT(flagMTset, flagDataset, flagRelaxation, flagMethod, flagLipidSet,listMTs_mpi,flagPruning);
        growthMTs = growthMTs_mpi;
        matchedPhenotypes = matchedPhenotypes_mpi;
        processedData = processedData_mpi;
        
        
    case {'lusk'}
        flagMTset='lusk';
        flagDataset=[];
        flagRelaxation=true;
        flagMethod='ACHR';
        flagLipidSet='lusk';
        [growthMTs_lusk, matchedPhenotypes_lusk, processedData_lusk] = simulateGrowthMT(flagMTset, flagDataset, flagRelaxation, flagMethod, flagLipidSet,listMTs_lusk,flagPruning);
        growthMTs = growthMTs_lusk;
        matchedPhenotypes = matchedPhenotypes_lusk;
        processedData = processedData_lusk; 
end


% select feasible models:
idxFeasible = ~cellfun(@isempty, growthMTs(:,2));
feasibleModels = growthMTs(idxFeasible,:);





% ###
% ##### STEP 2: Analyze the results and save into struct:
% ###

switch flagMethod
    case {'ACHR'}
        % load sampling model for WT:
        dirModels = dir(fullfile('..','SamplingResults','Models_sampling','modelSamplingWT.mat'));
        modelSamplingWT = load(fullfile(dirModels.folder,'modelSamplingWT.mat'));
        modelSamplingWT = modelSamplingWT.modelSamplingWT;
        
        % get flux samples for WT:
        dirSamples = dir(fullfile('..','SamplingResults','Random_samples','tempVariable.mat'));
        pathSampless = dirSamples.folder;
        
        % load samples and model for WT:
        rootSamplesFile = 'samplesWT_ACHR_';
        samplesWT = zeros(numel(modelSamplingWT.rxns),5000);
        count = 1;
        for i = 1:5
            nameSamplesFile_i = [rootSamplesFile, char(string(i)), '.mat'];
            pathSamples_i = fullfile(pathSampless, nameSamplesFile_i);
            samplesWT_i = load(pathSamples_i);
            samplesWT_i = samplesWT_i.samplesWT_i;
            samplesWT(1:size(samplesWT_i,1), count:count+(size(samplesWT_i,2)-1)) = samplesWT_i;
            count = count+(size(samplesWT_i,2));
        end    
        
        
        % get fluxes sampled for mutant lines:
        fluxesKOs = '';
        dirSampledFluxes = dir(fullfile('..','SamplingResults','SampledFluxes','tempVariable.mat'));
        pathSampledFluxes = dirSampledFluxes.folder;
        files_wRelax = cell(size(growthMTs,1),1);
        rootFolder = strrep(pathSampledFluxes, '\SampledFluxes', '');
        
        for i = 2:size(feasibleModels,1)
            ko_i = feasibleModels{i,1};
            fileKO_i = ['ACHR_wRelax_',ko_i,'.txt'];
            pathKO_i = fullfile(pathSampledFluxes,fileKO_i);
            files_wRelax{i} = ko_i;
            
            % get list of reactions for WT and retrieve corresponding sampled fluxes:
            idxRxnsWT = contains(modelSamplingWT.grRules, ko_i);
            fluxesWT = samplesWT(idxRxnsWT,:);
            
            % get fluxes sampled for selected KO-rxn with relaxation of bounds:
            sampledFluxesKOi_w = getFluxes4KOrxns(ko_i, fluxesWT, 'w_relax', flagMethod, rootFolder);
            
            % save variable into struct:
            fluxesKOs.(ko_i).wRelaxation = sampledFluxesKOi_w;

            % export results as .txt file:
            writetable(sampledFluxesKOi_w, pathKO_i)

        end
        
end

path2save_wRelax = fullfile(pathSampledFluxes,'ACHR_wRelax.txt');
idxEmpty = cellfun(@isempty, files_wRelax);
files_wRelax(idxEmpty) = '';
writecell(files_wRelax, path2save_wRelax)





% ###
% ##### STEP 3: Compare the means between WT (Col-0) and MTs (T-DNA lines):
% ###
dirTtests = dir(fullfile('..','SamplingResults','Ttest','tempVariable.mat'));
pathStats = fullfile(dirTtests.folder, 'Results_statistic_tests.xlsx');

% classify mutants according to the GPR rule of the corresponding
% reactions: (i) the reaction is catalyzed by a unique enzyme, or (ii) the
% reaction is catalyzed by a set of isoenzymes
listUniqueEnzymes = matchedPhenotypes(strcmp(matchedPhenotypes(:,3), '0'),1);
listIsoenzymes = matchedPhenotypes(strcmp(matchedPhenotypes(:,3), '1'),1);


% classify the mutants corresponding to reactions catalyzed by isoenzymes
% according to differences in means using results from T-tests: 
% NOTE: For this the ratio between mean WT/MT is calculated and the
% reactions will be categorized as follows:
% (i) Full rescue: reactions in which mean-WT == mean-MT
% (ii) Partial rescue: reactions in which the mean-MT < mean-WT, but mean-MT ~= 0
% (iii) Full-KO: reactions for which mean-MT == 0.

pathClassifiedMTs = fullfile(dirTtests.folder, 'Results_classification_mutants.xlsx');

% carry out T-test for dataset with relaxation, to determine if the
% means are different from zero (0) and if mean_WT is different from mean_MT:
if ~isempty(cellfun(@isempty, files_wRelax))
    listMTs_wRelax = files_wRelax;
    [meanMTs_wR, oneSampleTtest_wR, twoSampleTtest_wR] = perform_Ttest(listMTs_wRelax, 'w_relax',pathSampledFluxes);
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
        
    % add phenotypes info:
    [~,idxPhenotype] = ismember(mutantsClasses_wR(:,1),matchedPhenotypes(:,1));
    idxPhenotype(1) = 1;
    idxPhenotype(idxPhenotype==0) = '';
    mutantsClasses_wR = [mutantsClasses_wR, matchedPhenotypes(idxPhenotype,:)];
    
    writecell(mutantsClasses_wR, pathClassifiedMTs, 'Sheet', 'mutantsClasses_wR')
end



end