function model = computefluxsumMT(processedData, averageFC_KOi, sdFC_KOi, referenceFluxes_wSFR, flagFieldType, solFVA)
%##########################################################################
%
% This function calculates the flux sum for the Arabidopsis mutants, by
% using the flux sum of reactions of the wild-type (v_i^Col0) as reference
% and the ratio of relative abundances measured on the MT and WT (m =
% RA-MT/RA-WT). In addition, the lower and upper bounds for the reactions
% are relaxed by considering the standard deviation (SD):
% For this the following formula will be applied:
%
%  v_i^MT <= v_i^Col0 * (m + SD)
%  v_i^MT >= v_i^Col0 * (m - SD)
% 
% flagRelaxation - (true) Use sd to calculate lower and upper bounds
%                  (false) Default. Do not relax upper and lower bounds
%
% flagFieldType - (0) Default. The incoming fluxes are stored in nested struct fields
%                 (1) The incoming fluxes are not stored as nested fields
%
%##########################################################################

if nargin < 5
    flagFieldType = 0;
end

% ###
% ##### 1). Create pseudo-metabolites and -reactions for each of the lipid pools:
% ###
model = processedData.modelIrrev;
model.subSystems = cellstr(string(model.subSystems));
listLipidPools = fieldnames(referenceFluxes_wSFR);

for i = 1:numel(listLipidPools)
    lipidPool_i = listLipidPools{i};
    
    
    % (i). get FSR for the selected pool:
    fsr_i = referenceFluxes_wSFR.(lipidPool_i).fsr;
    
    
    % (ii). get indexes for all reactions contributing to the pool:
    switch flagFieldType
        case 0
            speciesPool_i = fieldnames(referenceFluxes_wSFR.(lipidPool_i));
            idxRxnsPool_i = zeros(100,1);
            count = 1;
            for j = 2:(numel(speciesPool_i)-1)
                idxsMet_j = referenceFluxes_wSFR.(lipidPool_i).(speciesPool_i{j}).idxIncomingRxns;
                idxRxnsPool_i(count:count+(numel(idxsMet_j)-1)) = idxsMet_j;
                count = count+numel(idxsMet_j);
            end
            idxRxnsPool_i(idxRxnsPool_i==0) = '';
        
        case 1
            idxRxnsPool_i = referenceFluxes_wSFR.(lipidPool_i).idxIncomingRxns;
    end
    idxRxnsPool_i = unique(idxRxnsPool_i);
    
    % (iii)add new pseudo-reaction:
    model.S(size(model.S,1),end+1) = 0;
    model.rxns(end+1)           = {['pseudoRxn_',lipidPool_i]};
    model.c(end+1)              = 0;
    model.rules(end+1)          = {''};
    model.grRules(end+1)        = {''};
    model.rxnGeneMat(end+1,:)   = 0;
    model.rxnConfidenceScores(end+1) = 0;
    model.rxnNames(end+1)       = {'pseudo-reaction'};
    model.rxnNotes(end+1)       = {''};
    model.rxnECNumbers(end+1)   = {''};
    model.rxnReferences(end+1)  = {''};
    model.rxnKEGGID(end+1)      = {''};
    model.subSystems(end+1)     = {''};
    model.match(end+1,:)        = 0;
    model.C                     = sparse(0,numel(model.rxns));
    
    % (iv)add new pseudo-metabolite:
    model.S(end+1,idxRxnsPool_i) = 1;
    model.S(end,end)             = -1;
    model.b(end+1)               = 0;
    model.mets(end+1)            = {['pseudoMet_',num2str(i)]};
    model.metNames(end+1)        = {'pseudo-metabolite'};
    model.csense(end+1)          = 'E';
    model.metCharges(end+1)      = 0;
    model.metFormulas(end+1)     = {''};
    model.metInChIString(end+1)  = {''};
    model.metKEGGID(end+1)       = {''};
    model.metChEBIID(end+1)      = {''};
    model.metPubChemID(end+1)    = {''};
    model.metLIPIDMAPSID(end+1)  = {''};
    
    
    % (v) constrain the LB and UB of the pseudo-reaction:
    fcPool_i = str2double(string(averageFC_KOi(i+1,2)));
    sdPool_i = str2double(string(sdFC_KOi(i+1,2)));
    
    ubPool_i = fsr_i*(fcPool_i+sdPool_i);
    model.ub(end+1) = ubPool_i;
    
    if nargin < 6 % compute lb using fsr and relaxing with sd
        lbPool_i = fsr_i*(fcPool_i-sdPool_i);
        
    else % relax the lb using the minimum value of the feasible space computed via FVA
        minFluxPool_i = solFVA.minFlux(idxRxnsPool_i); % get minFlux for reactions
        lbPool_i = min(minFluxPool_i);
    end
    model.lb(end+1) = lbPool_i;
end

end

