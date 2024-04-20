function modelMT = createModelMT(processedData, idxKO, referenceFluxes_wSFR)
%##########################################################################
%
% This function is used to create the model for Arabidopsis mutants (MT) by
% constraining lipid-related reactions with the flux_sum_MT. The latter is
% calculated by adjusting a reference flux sum obtained for Arabidopsis
% wild-type (WT) with the ratios of relative abundance of lipids among 
% MT/WT. 
%
% INPUT:
%   processedProfiles - ratios of relative abundances (RA) among MT/WT
%   idxKO             - index of mutant for which lipid ratios
%                       must be retrieved to constrain the model
%
% OUTPUT: modelMT     - COBRA model constrained using the flux sum
%                       calculated for the selected Arabidopsis mutant
%
%##########################################################################

processedProfiles = processedData.processedProfiles;

% get FC calculated using average relative abundance:
averageFC_KOi = processedProfiles.FCaverages(:,idxKO);

% get SD calculated for the fold-change of the biological replicates:
sdFC_KOi = processedProfiles.SD_FC(:,idxKO);


% compute flux sum for MT and constrain boounds of lipid-related rxns:
modelMT = computefluxsumMT(processedData, averageFC_KOi, sdFC_KOi, referenceFluxes_wSFR);

