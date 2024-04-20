function [meanMTs, oneSampleTtest, twoSampleTtest] = perform_Ttest(listMTs, flagSet, pathSampledFluxes)
%##########################################################################
%
%
% flagSet   - ('wo_relax') without relaxation of lower and upper bounds using lipids sd
%           - ('w_relax') with relaxation of lower and upper bounds using lipids sd
% 
% The One-sample T-test is implemented to test if the means of random
% fluxes are different from zero. Hence, Ho: µ=0; Hi: µ?0.
%    
% The Two-sample T-test is implemented to test that means of random fluxes
% between WT and MT are different from each other. Hence, Ho: µ1=µ2; Hi: µ1?µ2
%
% If the value returned for h is h = 1 indicates that the T-test rejects the
% null hypothesis at the default significance level of 5%.
%
%##########################################################################


switch flagSet
    case {'wo_relax'}
        prefixFile = 'ACHR_woRelax_';

    case {'w_relax'}
        prefixFile = 'ACHR_wRelax_';
end

listMT_IDs = cell(size(listMTs,1),1);
listMT_IDs{1} = 'mutantID';
meanMTs = zeros(size(listMTs,1),4);% 'meanWT' 'varianceWT' 'meanMT' 'varianceMT'
twoSampleTtest = zeros(numel(listMTs),8);% 'h' 'p' 'ci_lb' 'ci_ub' 'tstat' 'df' 'sd_lb' 'sd_ub'
labelsTwoSampleTest = cell(size(twoSampleTtest,1),1);
oneSampleTtest = zeros(numel(listMTs)*2,9);% 'h' 'p' 'ci_lb' 'ci_ub' 'tstat' 'df' 'sd_lb' 'sd_ub' 'No. Observ'
labelsRowsTtest = cell(numel(listMTs)*2,2);
labelsColsTtest = {'locusID','dataSet', 'h', 'p', 'ci_lb', 'ci_ub', 'tstat', 'df', 'sd_lb', 'sd_ub', 'No. Observ'};
count = 1;

% calculate means for random fluxes of WT and MTs:
for i = 1:size(listMTs,1)
    % load file with sampled fluxes:
    ko_i = listMTs{i};
    nameFile_MTi = [prefixFile, ko_i, '.txt'];
    samplesKO_i = readcell(fullfile(pathSampledFluxes,nameFile_MTi));
    
    % split samples into WT and MT:
    samplesWT_i  = samplesKO_i(strcmp(samplesKO_i(:,1), 'WT'),2);
    samplesWT_i = str2double(string(samplesWT_i));
    
    samplesMT_i  = samplesKO_i(strcmp(samplesKO_i(:,1), 'MT'),2);
    samplesMT_i = str2double(string(samplesMT_i));
    
    
    % calculate media and variance for each data set and save results:
    listMT_IDs{i} = ko_i;
    meanMTs(i,1) = mean(samplesWT_i);
    meanMTs(i,2) = var(samplesWT_i);
    meanMTs(i,3) = mean(samplesMT_i);
    meanMTs(i,4) = var(samplesMT_i);
    
    
    % Implement One-sample T-test to test that means are different from
    % zero: Ho: µ=0; Hi: µ?0
    labelsRowsTtest(count:count+1,1) = cellstr(string(ko_i));
    labelsRowsTtest(count:count+1,2) = {'WT';'MT'};
    
    [h1wt,p1wt,ci1wt,stats1wt] = ttest(samplesWT_i);
    oneSampleTtest(count,1:6) = [h1wt,p1wt,ci1wt(1),ci1wt(2),stats1wt.tstat,stats1wt.df];
    oneSampleTtest(count,9) = numel(samplesWT_i);
    if numel(stats1wt.sd) == 1
        oneSampleTtest(count,7) = [stats1wt.sd];
    else
        oneSampleTtest(count,7:8) = [stats1wt.sd(1),stats1wt.sd(2)];
    end
    
    [h1mt,p1mt,ci1mt,stats1mt] = ttest(samplesMT_i);
    oneSampleTtest(count+1,1:6) = [h1mt,p1mt,ci1mt(1),ci1mt(2),stats1mt.tstat,stats1mt.df];
    oneSampleTtest(count+1,9) =  numel(samplesMT_i);
    if numel(stats1mt.sd) == 1
        oneSampleTtest(count+1,7) = [stats1mt.sd];
    else
        oneSampleTtest(count+1,7:8) = [stats1mt.sd(1),stats1mt.sd(2)];
    end
    count = count+2;
    
    
    % Implement Two-sample T-test to test that means between WT and MT are
    % different from each other: Ho: µ1=µ2; Hi: µ1?µ2
    
    % check if number of observations is the same and act accordingly:
    if numel(samplesWT_i) == numel(samplesMT_i)
        % implement equal variance (independent) T-test:
        [h2,p2,ci2,stats2] = ttest2(samplesWT_i, samplesMT_i);
        labelsTwoSampleTest{i} = 'equal variance (independent)';
    else
        % check variances ratios:
        if meanMTs(i,2) > meanMTs(i,4)
            diffVariances = meanMTs(i,2)/meanMTs(i,4);
        else
            diffVariances = meanMTs(i,4)/meanMTs(i,2);
        end
        
        % to determine if variances are different, apply rule of thumb: if
        % the ratio of the larger variance to the smaller variance is less
        % than 4 then we can assume the variances are approximately equal
        if diffVariances < 4
            % implement equal variance (independent) T-test:
            [h2,p2,ci2,stats2] = ttest2(samplesWT_i, samplesMT_i);
            labelsTwoSampleTest{i} = 'equal variance (independent)';
            
        else
            % implement unequal variance (independent) T-test:
            [h2,p2,ci2,stats2] = ttest2(samplesWT_i, samplesMT_i, 'Vartype','unequal');
            labelsTwoSampleTest{i} = 'unequal variance (independent)';
        end
    end
    
    % add the results to the table:
    twoSampleTtest(i,1:6) = [h2,p2,ci2(1),ci2(2),stats2.tstat,stats2.df];
    
    if numel(stats2.sd) == 1
        twoSampleTtest(i,7) = [stats2.sd];
    else
        twoSampleTtest(i,7:8) = [stats2.sd(1),stats2.sd(2)];
    end
end

% add labels to statistics results:
meanMTs = cellstr(string(meanMTs));
meanMTs = [listMT_IDs, meanMTs];
meanMTs = [{'locusID','meanWT','varianceWT','meanMT','varianceMT'};meanMTs];

oneSampleTtest = [labelsRowsTtest, cellstr(string(oneSampleTtest))];
oneSampleTtest = [labelsColsTtest;oneSampleTtest];

twoSampleTtest = [listMTs, labelsTwoSampleTest, cellstr(string(twoSampleTtest))];
labelsColsTtest{2} = 'T-test type';
twoSampleTtest = [labelsColsTtest(1:numel(labelsColsTtest)-1);twoSampleTtest];


end

  