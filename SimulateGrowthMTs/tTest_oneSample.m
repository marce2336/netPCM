function oneSampleTtest = tTest_oneSample(oneSampleTtest,samplesWT_i,count)
%##########################################################################
%
% The One-sample T-test is implemented to test if the means of random
% fluxes are different from zero. Hence, Ho: µ=0; Hi: µ?0.
%    
% If the value returned for h is h = 1 indicates that the T-test rejects the
% null hypothesis at the default significance level of 5%.
%
%##########################################################################


[h1wt,p1wt,ci1wt,stats1wt] = ttest(samplesWT_i);
oneSampleTtest(count,1:6) = [h1wt,p1wt,ci1wt(1),ci1wt(2),stats1wt.tstat,stats1wt.df];
oneSampleTtest(count,9) = numel(samplesWT_i);
if numel(stats1wt.sd) == 1
    oneSampleTtest(count,7) = [stats1wt.sd];
else
    oneSampleTtest(count,7:8) = [stats1wt.sd(1),stats1wt.sd(2)];
end

end

