function [models, samples, resultsTtests, summaryFCs, summaryGeneExpression] = differentialFluxAnalysis
%##########################################################################
%
% This script is used to find out if the mutants for which absence of 
% phenotypic changes is predicted, may exhibit differences in other 
% pathways indirectly or not related to the metabolic steps affected by the
% mutation.
%
%##########################################################################

fullRescueMTs = {'AT1G06520outl', 'AT1G02390outl','AT3G14360', 'AT3G03540'};

% ###
% ##### Step 1: load models and samples
% ###
dirModels = dir(fullfile('..','SamplingResults','Models_sampling','tempVariable.mat'));
pathModels = dirModels.folder;
pathSamples = strrep(pathModels,'Models_sampling','Random_samples');

for i = 1:numel(fullRescueMTs)
    mt_i = fullRescueMTs{i};
    modelName = [mt_i, '_modelSampling_w_relax.mat'];
    modelMTi = load(fullfile(pathModels,modelName));
    modelMTi = modelMTi.modelSampling;
    models.(mt_i) = modelMTi;
    
    % load samples:
    newFlag = 'w_relax_';
    samplesMTi = zeros(numel(modelMTi.rxns),5000);
    count = 1;
    for j = 1:5
        nameSamplesFile_j = [mt_i, '_ACHR_', newFlag, char(string(j)), '.mat'];
        pathSamples_j = fullfile(pathSamples, nameSamplesFile_j);
        samplesMT_j = load(pathSamples_j);
        samplesMT_j = samplesMT_j.sSamples_j;
        samplesMTi(1:size(samplesMT_j,1), count:count+(size(samplesMT_j,2)-1)) = samplesMT_j;
        count = count+(size(samplesMT_j,2));
    end
    samples.(mt_i) = samplesMTi;
end




% ###
% ##### Step 2: compute means for each reaction and perform two-sample t-Test
% ###

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
    samplesWT(1:size(samplesWT_i,1), count:count+(size(samplesWT_i,2)-1)) = samplesWT_i;
    count = count+(size(samplesWT_i,2));
end


resultsTtests = '';

for i = 1:numel(fullRescueMTs)
    mt_i = fullRescueMTs{i};
    samplesMT_i = samples.(mt_i);
    
    modelMTi = models.(mt_i);
    rxnsMT_i = modelMTi.rxns;
    
    % remove from WT samples, the reactions that were eliminated during
    % model reduction of MT:
    rxnsWT = modelSamplingWT.rxns;
    [idxRxns2Keep, locb] = ismember(rxnsWT, rxnsMT_i);
    rxns2Keep = rxnsWT(idxRxns2Keep);
    prunedSamplesWT = samplesWT(idxRxns2Keep,:);
    prunedSubsystems = modelSamplingWT.subSystems(idxRxns2Keep);
    
    % pre-allocate space to save t-Test results:
    oneSampleTest_MTi = zeros(sum(idxRxns2Keep),10);
    %labelsCols = {'h', 'p', 'ci_lb', 'ci_ub', 'tstat', 'df', 'sd_lb', 'sd_ub', 'No. Observ', 'meanSamples'};
    twoSampleTest_MTi = zeros(sum(idxRxns2Keep),8);
    %labelsCols = {'h' 'p' 'ci_lb' 'ci_ub' 'tstat' 'df' 'sd_lb' 'sd_ub'};
    
    
    % perform one- and two-sample t-Test:
    for j = 1:numel(rxns2Keep)
        rxn_j = rxns2Keep{j};
        idxRxnj_MT = strcmp(modelMTi.rxns, rxn_j);
        samplesMTi_rxnj = samplesMT_i(idxRxnj_MT,:);
        
        % compute mean of fluxes:
        meanMTi_Rnxj = mean(samplesMTi_rxnj);
        oneSampleTest_MTi(j, end) = meanMTi_Rnxj;
    
        % carry one-sample t-Test for each reaction:
        oneSampleTest_MTi = tTest_oneSample(oneSampleTest_MTi,samplesMTi_rxnj,j);
        
        % carry out two-sample t-Test to find differences among means of WT and MT_i:
        samplesWT_rxnj = prunedSamplesWT(j,:);
        [twoSampleTest_MTi,~] = tTest_twoSample(samplesWT_rxnj,samplesMTi_rxnj,twoSampleTest_MTi,j);
    end
    resultsTtests.(mt_i).oneSampleTtest = oneSampleTest_MTi;
    resultsTtests.(mt_i).twoSampleTtest = twoSampleTest_MTi;
    resultsTtests.(mt_i).rxns = rxns2Keep;
    resultsTtests.(mt_i).subSystems = prunedSubsystems;
    resultsTtests.(mt_i).prunedSamplesWT = prunedSamplesWT;
    
    % save samples for MT:
    locb = locb(locb~=0);
    prunedSamplesMT = samplesMT_i(locb,:);
    resultsTtests.(mt_i).samplesMT = prunedSamplesMT;
end




% ###
% ##### Step 3: compute pathways enriched in differential fluxes
% ###

for i = 1:numel(fullRescueMTs)
    mt_i = fullRescueMTs{i};
    
    % find number of reactions with mean of fluxes equal to WT:
    twoTtestMT_i = resultsTtests.(mt_i).twoSampleTtest;
    idxRxnsSteadyFlux = twoTtestMT_i(:,1) == 0;
    subSystemsMT_i = resultsTtests.(mt_i).subSystems;
    subsystemsSteadyMT_i = subSystemsMT_i(idxRxnsSteadyFlux);
    uniqueSubsystems = unique(cellstr(string(subsystemsSteadyMT_i)));

    % count rxns with equal mean fluxes per subsystem:
    equalFluxes = zeros(numel(uniqueSubsystems),1);
     for j = 1:numel(uniqueSubsystems)
         idxHits = strcmp(subsystemsSteadyMT_i, uniqueSubsystems{j});
         equalFluxes(j) = sum(idxHits);
     end
     [equalFluxes, locb] = sort(equalFluxes, 'descend');
     equalFluxes = [uniqueSubsystems(locb),cellstr(string(equalFluxes))];
     equalFluxes = [{'subSystem','No. reactions'};equalFluxes];
     resultsTtests.(mt_i).equalFluxes = equalFluxes;
     
     
     % find number of reactions with differential fluxes compared to WT:
    subSystemsMT_i = resultsTtests.(mt_i).subSystems;
    rxnsMT_i = resultsTtests.(mt_i).rxns;
    subsystemsSteadyMT_i = subSystemsMT_i(idxRxnsSteadyFlux==0);
    rxnsMT_i = rxnsMT_i(idxRxnsSteadyFlux==0);
    uniqueSubsystems = unique(cellstr(string(subsystemsSteadyMT_i)));
    idxEmpty = cellfun(@isempty, uniqueSubsystems);
    uniqueSubsystems(idxEmpty) = '';
    listDiffRxns.uniqueSubsystems = uniqueSubsystems;
    
    % count rxns with differential fluxes per subsystem:
    differentialFluxes = zeros(numel(uniqueSubsystems),1);
     for j = 1:numel(uniqueSubsystems)
         idxHits = strcmp(subsystemsSteadyMT_i, uniqueSubsystems{j});
         differentialFluxes(j) = sum(idxHits);
         pseudoName_j = ['subsystem_', char(string(j))];
         listDiffRxns.(pseudoName_j).rxns = rxnsMT_i(idxHits);
     end
     [differentialFluxes, locb] = sort(differentialFluxes, 'descend');
     differentialFluxes = [uniqueSubsystems(locb),cellstr(string(differentialFluxes))];
     differentialFluxes = [{'subSystem','No. reactions'};differentialFluxes];
     list2remove = {'SLIME pseudo reaction';'transport';'export';'demand';'import';'biomass reaction';'maintenance';'im-/export';'Exchange'};
     idx2remove = ismember(differentialFluxes(:,1), list2remove);
     differentialFluxes(idx2remove,:) = '';
     resultsTtests.(mt_i).differentialFluxes = differentialFluxes;
     resultsTtests.(mt_i).listDiffRxns = listDiffRxns;
     
     % filter reactions with differential fluxes:
     differentialRxns = resultsTtests.(mt_i).rxns(idxRxnsSteadyFlux==0);
     differentialSubsystems = subsystemsSteadyMT_i;
     differentialSamples_MT = resultsTtests.(mt_i).samplesMT(idxRxnsSteadyFlux==0,:);
     differentialSamples_WT = resultsTtests.(mt_i).prunedSamplesWT(idxRxnsSteadyFlux==0,:);
     
     resultsTtests.(mt_i).differentialRxns = differentialRxns;
     resultsTtests.(mt_i).differentialSubsystems = differentialSubsystems;
     resultsTtests.(mt_i).differentialSamples_MT = differentialSamples_MT;
     resultsTtests.(mt_i).differentialSamples_WT = differentialSamples_WT;
end




% ###
% ##### Step 4: analyze magnitude of fold changes for differential reactions
% ###

for i = 1:numel(fullRescueMTs)
    mt_i = fullRescueMTs{i};
    differentialFluxes_mti = resultsTtests.(mt_i).differentialFluxes;
    maxRxns_mti = max(str2double(string(differentialFluxes_mti(2:end,2))));
    fc_mti = zeros(size(differentialFluxes_mti,1),maxRxns_mti);
    fcAlt_mti = fc_mti;
    groupsFC_mti = zeros(size(differentialFluxes_mti,1),5);
    barplotFC_mti = zeros((size(differentialFluxes_mti,1)-1)*4,1);
    rowNamesBarplot_mti = cell(numel(barplotFC_mti),1);
    count = 1;
    for j = 2:size(differentialFluxes_mti,1)
        % get subset of samples for selected subsystem:
        subsystem_j = differentialFluxes_mti{j,1};
        idxSubsystem_j = strcmp(resultsTtests.(mt_i).differentialSubsystems, subsystem_j);
        samplesSubj_wt = resultsTtests.(mt_i).differentialSamples_WT(idxSubsystem_j,:);
        samplesSubj_mt = resultsTtests.(mt_i).differentialSamples_MT(idxSubsystem_j,:);
    
        % compute fc: Ratio mean-flux-WT/-MT
        meanSubj_wt = mean(samplesSubj_wt,2);
        meanSubj_mt = mean(samplesSubj_mt,2);
        fcSubj = meanSubj_wt./meanSubj_mt;
        fcSubj_alt = meanSubj_mt./meanSubj_wt;
        
        fc_mti(j,1:numel(fcSubj)) = fcSubj';
        fcAlt_mti(j,1:numel(fcSubj)) = fcSubj_alt';
        
        % count reactions per subsystem within range provided:
        rowNamesBarplot_mti(count:count+3) = cellstr(string(subsystem_j));
        prunedFCsubj = fcSubj(fcSubj~=0);
        idxCount = prunedFCsubj < 0.05;
        groupsFC_mti(j,1) = sum(idxCount);
        barplotFC_mti(count) = sum(idxCount);
        
        idxCount = (prunedFCsubj < 0.5).*(prunedFCsubj > 0.05);
        groupsFC_mti(j,2) = sum(idxCount);
        barplotFC_mti(count+1) = sum(idxCount);
        
        idxCount = (prunedFCsubj < 50).*(prunedFCsubj > 2);
        groupsFC_mti(j,3) = sum(idxCount);
        barplotFC_mti(count+2) = sum(idxCount);
        
        idxCount = prunedFCsubj > 50;
        groupsFC_mti(j,4) = sum(idxCount);
        barplotFC_mti(count+3) = sum(idxCount);
        count = count+4;
        groupsFC_mti(j,5) = numel(prunedFCsubj);
    end
    
    fc_mti = [differentialFluxes_mti(:,1), cellstr(string(fc_mti))];
    fcAlt_mti = [differentialFluxes_mti(:,1), cellstr(string(fcAlt_mti))];
    resultsTtests.(mt_i).fc = fc_mti;
    resultsTtests.(mt_i).fc_alternative = fcAlt_mti;
    resultsTtests.(mt_i).groupsFC = groupsFC_mti;

    % convert categorization of FC into table and save as struct:
    filterNames = cell2table(repmat({'< 0.05';'0.05< <0.5';'2< <50';'>50'},numel(barplotFC_mti)/4,1));
    filterNames.Properties.VariableNames = "range";
    rowNamesBarplot_mti = cell2table(rowNamesBarplot_mti);
    rowNamesBarplot_mti.Properties.VariableNames = "subsystem";
    barplotFC_mti = array2table(barplotFC_mti);
    barplotFC_mti.Properties.VariableNames = "count";
    barplotFC_mti = [rowNamesBarplot_mti, barplotFC_mti, filterNames];
    resultsTtests.(mt_i).barplotFC = barplotFC_mti;
end




% ###
% ##### Step 5: group subsystems into metabolite categories
% ###

path2save = fullfile(strrep(pathModels,'Models_sampling','Ttest'),'barplots_differential_fluxes.xlsx');
pathCategories = fullfile('InputDataSimulateGrowth','List_subsystem_categories.xlsx');
listCategories = cellstr(string(readcell(pathCategories, 'Sheet', 'list_categories'))); %list_abbreviations
uniqueCategories = unique(listCategories(2:end,3));
listRanges = {'< 0.05';'0.05< <0.5';'2< <50';'>50'};

% create empty table:
a = cell2table({'mets'});
a.Properties.VariableNames = "metabolite_class";
b = array2table(0);
b.Properties.VariableNames = "count";
c = cell2table({'range'});
c.Properties.VariableNames = "range";
d = array2table([0,0]);
d.Properties.VariableNames = ["total","proportion"];
e = cell2table({'locus'});
e.Properties.VariableNames = "locus";
mergedBarplot = [a,b,c,d,e];



for i = 1:numel(fullRescueMTs)
    mt_i = fullRescueMTs{i};
    barPlot_mti = resultsTtests.(mt_i).barplotFC;
    subsystems_mti = table2cell(barPlot_mti(:,1));
    ranges_mti = table2cell(barPlot_mti(:,3));
    summaryBarplot_mti = cell(numel(uniqueCategories)*4,5);
    count = 1;
    
    for j = 1:numel(uniqueCategories)
        metCategory_j = uniqueCategories{j};
        idxCategory_j = strcmp(listCategories(:,3), metCategory_j);
        subsCategory_j = listCategories(idxCategory_j,1);
        
        for k = 1:numel(listRanges)
            range_k = listRanges{k};
            idxMetsRange_k = ismember(subsystems_mti, subsCategory_j).*strcmp(ranges_mti, range_k);
            sumRange_k = sum(table2array(barPlot_mti(idxMetsRange_k==1,2)));
            summaryBarplot_mti{count,1} = metCategory_j;
            summaryBarplot_mti{count,2} = cellstr(string(sumRange_k));
            summaryBarplot_mti{count,3} = unique(ranges_mti(idxMetsRange_k==1));
            
            % compute total reactions within metabolite category:
            idxTotalCategori_k = ismember(subsystems_mti, subsCategory_j);
            totalRxnsCategory = sum(table2array(barPlot_mti(idxTotalCategori_k,2)));
            summaryBarplot_mti{count,4} = cellstr(string(totalRxnsCategory));
            
            % compute proportion of reactions whithin given range:
            proportionMetsRange_k = sumRange_k/totalRxnsCategory;
            summaryBarplot_mti{count,5} = cellstr(string(proportionMetsRange_k));
            count = count+1;
            
        end

    end
    
    % convert summary bar plot into table and save as struct:
    metClass_mti = cell2table(summaryBarplot_mti(:,1));
    metClass_mti.Properties.VariableNames = "metabolite_class";
    count_mti = array2table(str2double(string(summaryBarplot_mti(:,2))));
    count_mti.Properties.VariableNames = "count";
    range_mti = cell2table(summaryBarplot_mti(:,3));
    range_mti.Properties.VariableNames = "range";
    prop_mti = array2table(str2double(string(summaryBarplot_mti(:,4:5))));
    prop_mti.Properties.VariableNames = ["total","proportion"];
    tableBarplot_mti = [metClass_mti,count_mti,range_mti,prop_mti];
    resultsTtests.(mt_i).summaryBarplot = tableBarplot_mti;
    writetable(tableBarplot_mti, path2save, 'Sheet', mt_i)
    
    % create table with results for all T-DNA lines:
    locus_mti = cellstr(string(repmat(mt_i, size(tableBarplot_mti,1),1)));
    locus_mti = strrep(locus_mti, 'outl','');
    locus_mti = cell2table(locus_mti);
    locus_mti.Properties.VariableNames = "locus";
    tableBarplot_mti = [tableBarplot_mti,locus_mti];
    mergedBarplot = [mergedBarplot;tableBarplot_mti];
    
end

% adjust contents of 'range' column to improve visualization:
mergedBarplot(1,:) = [];
ranges = table2cell(mergedBarplot(:,3));
ranges = strrep(ranges, '< 0.05', 'a_< 0.05');
ranges = strrep(ranges, '0.05< <0.5', 'b_0.05< <0.5');
ranges = strrep(ranges, '2< <50', 'c_2< <50');
ranges = strrep(ranges, '>50', 'd_>50');
mergedBarplot(:,3) = cell2table(ranges);


% save summary metabolite class table in selected location:
writetable(mergedBarplot, path2save, 'Sheet', 'merged')




%--------------------------------------------------------------------------
% ###
% ##### Step 6: group subsystems of lipids network into metabolite categories
% ###
% this step is implemented to improve visualization
listLipidCategories = cellstr(string(readcell(pathCategories, 'Sheet', 'lipids_categories'))); %list_abbreviations
uniqueCategories = unique(listLipidCategories(2:end,2));
listRanges = {'< 0.05';'0.05< <0.5';'2< <50';'>50'};

% create empty table:
mergedLipidBarplot = [a,b,c,d,e];


for i = 1:numel(fullRescueMTs)
    mt_i = fullRescueMTs{i};
    barPlot_mti = resultsTtests.(mt_i).barplotFC;
    summaryLipidBarplot_mti = cell(numel(uniqueCategories)*4,5);
    subsystems_mti = table2cell(barPlot_mti(:,1));
    ranges_mti = table2cell(barPlot_mti(:,3));
    count = 1;
    
    for j = 1:numel(uniqueCategories)
        metCategory_j = uniqueCategories{j};
        idxCategory_j = strcmp(listLipidCategories(:,2), metCategory_j);
        subsCategory_j = listLipidCategories(idxCategory_j,1);
        
        for k = 1:numel(listRanges)
            range_k = listRanges{k};
            idxMetsRange_k = ismember(subsystems_mti, subsCategory_j).*strcmp(ranges_mti, range_k);
            sumRange_k = sum(table2array(barPlot_mti(idxMetsRange_k==1,2)));
            summaryLipidBarplot_mti{count,1} = metCategory_j;
            summaryLipidBarplot_mti{count,2} = cellstr(string(sumRange_k));
            summaryLipidBarplot_mti{count,3} = unique(ranges_mti(idxMetsRange_k==1));
            
            % compute total reactions whithin metabolite category:
            idxTotalCategori_k = ismember(subsystems_mti, subsCategory_j);
            totalRxnsCategory = sum(table2array(barPlot_mti(idxTotalCategori_k,2)));
            summaryLipidBarplot_mti{count,4} = cellstr(string(totalRxnsCategory));
            
            % compute proportion of reactions whithin given range:
            proportionMetsRange_k = sumRange_k/totalRxnsCategory;
            summaryLipidBarplot_mti{count,5} = cellstr(string(proportionMetsRange_k));
            
            count = count+1;
        end

    end
    
    % convert summary barplot into table and save as struct:
    metClass_mti = cell2table(summaryLipidBarplot_mti(:,1));
    metClass_mti.Properties.VariableNames = "metabolite_class";
    count_mti = array2table(str2double(string(summaryLipidBarplot_mti(:,2))));
    count_mti.Properties.VariableNames = "count";
    range_mti = cell2table(summaryLipidBarplot_mti(:,3));
    range_mti.Properties.VariableNames = "range";
    prop_mti = array2table(str2double(string(summaryLipidBarplot_mti(:,4:5))));
    prop_mti.Properties.VariableNames = ["total","proportion"];
    tableLipidBarplot_mti = [metClass_mti,count_mti,range_mti,prop_mti];
    resultsTtests.(mt_i).summaryBarplot = tableLipidBarplot_mti;
    nameSheet = ['lipid_',mt_i];
    writetable(tableLipidBarplot_mti, path2save, 'Sheet', nameSheet)
    
    % create table with results for all T-DNA lines:
    locus_mti = cellstr(string(repmat(mt_i, size(tableLipidBarplot_mti,1),1)));
    locus_mti = strrep(locus_mti, 'outl','');
    locus_mti = cell2table(locus_mti);
    locus_mti.Properties.VariableNames = "locus";
    tableLipidBarplot_mti = [tableLipidBarplot_mti,locus_mti];
    mergedLipidBarplot = [mergedLipidBarplot;tableLipidBarplot_mti];
      
end

% save summary lipid table in selected location:
mergedLipidBarplot(1,:) = [];
writetable(mergedLipidBarplot, path2save, 'Sheet', 'full_lipids')

% adjust contents of 'range' column to improve visualization:
ranges = table2cell(mergedLipidBarplot(:,3));
ranges = strrep(ranges, '< 0.05', 'a_< 0.05');
ranges = strrep(ranges, '0.05< <0.5', 'b_0.05< <0.5');
ranges = strrep(ranges, '2< <50', 'c_2< <50');
ranges = strrep(ranges, '>50', 'd_>50');
mergedLipidBarplot(:,3) = cell2table(ranges);


% eliminate uninteresting subsystems:
listSubsystems = table2cell(mergedLipidBarplot(:,1));
idxSubsystems = ismember(listSubsystems, {'lipoylated protein synthesis','oxylipin metabolism'});
mergedLipidBarplot(idxSubsystems,:) = [];
writetable(mergedLipidBarplot, path2save, 'Sheet', 'full_lipids_pruned')




%--------------------------------------------------------------------------
% ###
% ##### Step 7: take a look at the DEGs published by Bai et al. for GPAT1 mutant:
% ###
% Note: in the data published by Bai et al. (Bai Y, et al. A GPAT1 Mutation
% in Arabidopsis Enhances Plant Height but Impairs Seed Oil Biosynthesis.
% Int J Mol Sci. 2021 Jan 14;22(2):785. doi: 10.3390/ijms22020785. PMID:
% 33466786; PMCID: PMC7829857), the fold change (FC) for transcripts was
% computed between gpat1 and wild type (MT/WT). Values presented in the
% 'Supplemental Table S2' correspond to the mean log2 fold ratios.

pathDEGs = fullfile('InputDataSimulateGrowth','DEGs_Bai.xlsx');
baiDEGs = cellstr(string(readcell(pathDEGs, 'Sheet', 'pruned_DEGs')));
locusDEGs = baiDEGs(2:end, 2);

% get model for GPAT1 mutant (AT1G06520):
modelGPAT1 = models.AT1G06520outl;
idxModelDEGs = ismember(modelGPAT1.genes, locusDEGs);
hitsDEGs = modelGPAT1.genes(idxModelDEGs);

% get list of rxns (and subsystems) catalyzed by filtered gene products:
listDEGsRxns = '';
for i = 1:numel(hitsDEGs)
    locusDEG_i = hitsDEGs{i};
    idxRxns = contains(modelGPAT1.grRules, locusDEG_i);
    listDEGsRxns.(locusDEG_i).rxns = [modelGPAT1.rxns(idxRxns), modelGPAT1.subSystems(idxRxns)];
end

% get fold change for the reactions selected:
% Note: the FC among fluxes for Lusk data set used below was computed 
% previously (see Step 4, fcAlt_mti) as mean-flux-MT/WT. The fluxes were 
% obtained from files saved after running the 'processSamplingResults' function
fcGPAT1_lusk = resultsTtests.AT1G06520outl.fc_alternative;
listDiffRxns_gpat1 = resultsTtests.AT1G06520outl.listDiffRxns;

fieldDEGs = fieldnames(listDEGsRxns);
summaryFCs = cell(numel(fieldDEGs), 9);

for i = 1:numel(fieldDEGs)
    locus_i = fieldDEGs{i};
    summaryFCs{i,1} = locus_i;
    
    % get fc from Bai:
    idxLocus_i = strcmp(baiDEGs(:,2), locus_i);
    fcBai_locusi = baiDEGs(idxLocus_i, 4);
    summaryFCs{i,4} = fcBai_locusi;
    
    % get list rxns affected by mutation:
    rxnsMT_i = listDEGsRxns.(locus_i).rxns;
    summaryFCs{i,2} = join(rxnsMT_i(:,1), '|');
    summaryFCs{i,3} = join(rxnsMT_i(:,2), '|');
    uniqueSubsystems_mti = unique(rxnsMT_i(:,2));
    mergedFCs = zeros(size(rxnsMT_i,1),1);
    count = 1;
    
    for j = 1:numel(uniqueSubsystems_mti)
        uSub_j = uniqueSubsystems_mti{j};% get subsystem
        rxns_uSubj = rxnsMT_i(strcmp(rxnsMT_i(:,2),uSub_j),1);%list of rxns within subsystem
        idx_iSubj = find(strcmp(listDiffRxns_gpat1.uniqueSubsystems, uSub_j));
        nameField_uSubj = ['subsystem_', char(string(idx_iSubj))];
        rxnsLusk_uSubj = listDiffRxns_gpat1.(nameField_uSubj).rxns;  
        
        % find matching subsystems in Lusk fc table:
        idxHit = strcmp(fcGPAT1_lusk(:,1),uSub_j);
        
        % find rxn with differential fluxes and get FC value:
        matchingRxns = find(ismember(rxnsLusk_uSubj, rxns_uSubj));
        matchingRxns = matchingRxns+1;
        fcLusk_uSubj = fcGPAT1_lusk(idxHit,matchingRxns);
        mergedFCs(count:(count+numel(fcLusk_uSubj))-1) = str2double(string(fcLusk_uSubj));
        count = count+numel(fcLusk_uSubj);
    end
    
    % sort values of FC obtained:
    idxUp = mergedFCs > 1;
    summaryFCs{i,5} = cellstr(string(sum(idxUp)));
    summaryFCs{i,6} = cellstr(string(mean(mergedFCs(idxUp))));
    
    idxDown = mergedFCs < 1;
    summaryFCs{i,7} = cellstr(string(sum(idxDown)));
    summaryFCs{i,8} = cellstr(string(mean(mergedFCs(idxDown))));
    
end

summaryFCs(:,9) = cellstr(string(log2(str2double(string(summaryFCs(:,8))))));
summaryFCs(strcmp(summaryFCs(:,1), 'AT3G20160'),:) = '';

% convert summary into table and save as .xlsx file:
path2save = fullfile(strrep(pathModels,'Models_sampling','Ttest'),'Validated_TDNAs_full_rescue.xlsx');
locus = cell2table(summaryFCs(:,1:3));
locus.Properties.VariableNames = ["locus","rxns_ID","subSystems"];
log2s = array2table(str2double(string(summaryFCs(:,4:9))));
log2s.Properties.VariableNames = ["log2_FC_Bai","No_Rxns_Lusk_Up","mean_Rxns_Lusk_Up","No_Rxns_Lusk_Down","mean_Rxns_Lusk_Down","log2_mean_Rxns_Lusk_down"];
summaryTableFCs = [locus,log2s];
writetable(summaryTableFCs, path2save, 'Sheet','validated_locus_full_rescue')


% collapse list of FCs:
prefix = cellstr(string([11:1:32]));
prefix = strcat(prefix,'_');
prefix = [prefix,prefix];
prefix = cell2table(strcat(prefix',[summaryFCs(:,1);summaryFCs(:,1)]));
prefix.Properties.VariableNames = "locus";
trFCs = summaryTableFCs(:,4);
trFCs.Properties.VariableNames = "log2_FC";
fluxFCs = summaryTableFCs(:,9);
fluxFCs.Properties.VariableNames = "log2_FC";
mergedFCs = [trFCs;fluxFCs];
dataType = cell2table([repmat({'transcripts'},numel(trFCs),1); repmat({'fluxes'},numel(trFCs),1)]);
dataType.Properties.VariableNames = "dataType";
collapsedFCs = [prefix,dataType,mergedFCs];
writetable(collapsedFCs, path2save, 'Sheet','collapsed_validated_locus')




% ###
% ##### Step 8: get data of gene expression in leaves:
% ###
listLocus = cellstr(string(summaryFCs(2:end,1)));

listTissues = {'Leaf blade, intermediate 1','Leaf blade, intermediate 2','Leaf blade of the mature leaf',...
    'Leaf blade of the young leaf','Whole mature leaf','Leaf vein, intermediate 2','Vein of the mature leaf',...
    'Vein of the senescent leaf'}; 

pathExGs = fullfile('InputDataSimulateGrowth','DEGs_TRAVA');

summaryGeneExpression = array2table(zeros(numel(listLocus), numel(listTissues)));

for i = 1:numel(listLocus)
    %load file with expression data for locus:
    locus_i = listLocus{i};
    nameFile_i = [locus_i,'.xls'];
    exDataLocus_i = readcell(fullfile(pathExGs,nameFile_i));
    
    % get gene expression values for selected tissues:
    idxTissues = ismember(exDataLocus_i(:,1), listTissues);
    filteredTrs = str2double(string(exDataLocus_i(idxTissues,2)));
    summaryGeneExpression(i,:) = array2table(filteredTrs');
end

listLocus = cell2table(listLocus);
listLocus.Properties.VariableNames = "locus";
summaryGeneExpression.Properties.VariableNames = string(listTissues);
summaryGeneExpression = [listLocus,summaryGeneExpression];
writetable(summaryGeneExpression, path2save, 'Sheet','summary_gene_expression')


end

