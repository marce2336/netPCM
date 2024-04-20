function [twoSampleTtest,labelsTwoSampleTest] = tTest_twoSample(samplesWT_i,samplesMT_i,twoSampleTtest,i)
%##########################################################################
%
% % The Two-sample T-test is implemented to test that means of random fluxes
% between WT and MT are different from each other. Hence, Ho: µ1=µ2; Hi: µ1?µ2
%
% If the value returned for h is h = 1 indicates that the T-test rejects the
% null hypothesis at the default significance level of 5%.
%
%##########################################################################

% check if number of observations is the same and act accordingly:

if numel(samplesWT_i) == numel(samplesMT_i)
    % implement equal variance (independent) T-test:
    [h2,p2,ci2,stats2] = ttest2(samplesWT_i, samplesMT_i);
    labelsTwoSampleTest{i} = 'equal variance (independent)';
else
    % check variances ratios:
    varWT_i = var(samplesWT_i);
    varMT_i = var(samplesMT_i);
    
    if varWT_i > varMT_i
        diffVariances = varWT_i/varMT_i;
    else
        diffVariances = varMT_i/varWT_i;
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