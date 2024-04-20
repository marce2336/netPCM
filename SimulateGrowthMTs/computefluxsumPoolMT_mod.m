function model = computefluxsumPoolMT_mod(processedData, referenceFluxes, flagFieldType, solFVA, flagUseRange)
%##########################################################################
%
% This function calculates the flux sum for the Arabidopsis mutants, by
% using the flux sum of reactions of the wild-type (v_i^Col0) as reference
% and the ratio of relative abundances measured on the MT and WT (m =
% RA-MT/RA-WT). In addition, the lower and upper bounds for the reactions
% are relaxed by considering the standard deviation (SD).
% This script is a variation of the 'computeFluxSumMT', where all the
% species which are provided are considered as a unique pool. Hence, a
% single pseudo-metabolite/reaction are created.
%
% For this the following formula will be applied:
%
%  v_i^MT <= v_i^Col0 * (m + SD)
%  v_i^MT >= v_i^Col0 * (m - SD)
% 
% The relaxation of the bounds is carried out by running the function: 
% 'model = constrainRxnBounds(model, getFCspecies_i, idxIncomingRxns_ij, incomingFluxes_ij)'
%
% flagFieldType - (0) Default. The incoming fluxes are stored in nested struct fields
%                 (1) The incoming fluxes are not stored as nested fields
%
% solFVA        - Struct with flux vectors computed via FVA
%
% flagUseRange  - (0) Default. Bounds are relaxed using FSR computed via pFBA
%                 (1) Bounds are relaxed using min/max values from sampling
%
%##########################################################################

if nargin < 3
    flagFieldType = 0;
end

% ###
% ##### 1). Create pseudo-metabolites and -reactions for each of the lipid pools
% ###
model = processedData.model;
model.subSystems = cellstr(string(model.subSystems));
listLipidPools = fieldnames(referenceFluxes);


% ###
% ##### 2). Exclude constraints for selected pools
% ###
idxExcluded = ismember(listLipidPools, processedData.exclude);
listLipidPools(idxExcluded) = '';


% ###
% ##### 3). Get indexes for incoming reactions and add constraints
% ###
mt_i = processedData.mt_i;
for i = 1:numel(listLipidPools)
    lipidPool_i = listLipidPools{i};
    idxRxnsPool_i = zeros(300,1);
    count = 1;
    
    % get SD and FC values for pool_i:
    sdPool = referenceFluxes.(lipidPool_i).sd_mti;
    fcPool = referenceFluxes.(lipidPool_i).averageFC;
    
    % get indexes for all reactions contributing to the pool:
    
    % (i). get FSR for the selected pool:
    fsr_i = referenceFluxes.(lipidPool_i).fsr;
    
    % (ii). get indexes for all reactions contributing to the pool:
    speciesPool_i = fieldnames(referenceFluxes.(lipidPool_i));
    
    switch flagFieldType
        case 0
            for j = 1:(numel(speciesPool_i)-3)
                idxsMet_j = referenceFluxes.(lipidPool_i).(speciesPool_i{j}).idxIncomingRxns;
                idxRxnsPool_i(count:count+(numel(idxsMet_j)-1)) = idxsMet_j;
                count = count+numel(idxsMet_j);
            end


        case 1
            idxsMet_j = referenceFluxes.(lipidPool_i).idxIncomingRxns;
            idxRxnsPool_i(count:count+(numel(idxsMet_j)-1)) = idxsMet_j;
            count = count+numel(idxsMet_j);
    end
    idxRxnsPool_i(idxRxnsPool_i==0) = '';
    idxRxnsPool_i = unique(idxRxnsPool_i);
    
    % (iii). check if the pool should be split:
    switch mt_i
    % 'AT3G44830': the TG pool on the products side must be split into TGs produced via PDAT and DGAT enzymes
    % 'AT5G66450': the TG pool on the pdcts side is split into DGs pduced via PAP and SFR enzymes.
        case {'AT3G44830','AT5G66450','AT5G03080','AT5G57190','AT3G51520','AT2G19450','AT5G42870','AT1G15080','AT4G25970','AT4G22340'} 
            % Note: The locus 'AT4G00550','AT1G74320' were eliminated.
            % Spliting constraints make the models infeasible!
            % get name of rxn classes contributing to pool: ,
            nameRxnsPool_i = model.rxns(idxRxnsPool_i);
            rxnSuffix = processedData.rxnSuffix;
            poolSuffix = processedData.poolSuffix;
            poolNames = unique(extractBefore(nameRxnsPool_i, rxnSuffix));
            
            % separate indexes of rxns according to class:
            idxClasses = zeros(numel(idxRxnsPool_i),numel(poolNames));
            for j = 1:numel(poolNames)
                idxPool_j = contains(nameRxnsPool_i, poolNames{j});
                idxClasses(1:sum(idxPool_j),j) = idxRxnsPool_i(idxPool_j);
                poolNames{j} = [poolNames{j}, poolSuffix];
            end
            
        case {'AT2G39290'}
            % get name of rxn classes contributing to pool:
            nameRxnsPool_i = model.rxns(idxRxnsPool_i);
            nameRxnsPool_i = strrep(nameRxnsPool_i, 'PGPS_h', 'PGPSh_h');
            nameRxnsPool_i = strrep(nameRxnsPool_i, 'PGPS_m1', 'PGPSm_m1');
            nameRxnsPool_i = strrep(nameRxnsPool_i, 'PGPS_m2', 'PGPSm_m2');
            
            rxnSuffix = processedData.rxnSuffix;
            poolSuffix = processedData.poolSuffix;
            poolNames = unique(extractBefore(nameRxnsPool_i, rxnSuffix));
            
            % separate indexes of rxns according to class:
            idxClasses = zeros(numel(idxRxnsPool_i),2);
            for j = 1:numel(poolNames)
                idxPool_j = contains(nameRxnsPool_i, poolNames{j});
                idxClasses(1:sum(idxPool_j),j) = idxRxnsPool_i(idxPool_j);
                poolNames{j} = [poolNames{j}, poolSuffix];
            end
            
        case {'AT4G23850'}
            % get name of rxn classes contributing to pool:
            nameRxnsPool_i = model.rxns(idxRxnsPool_i);
            rxnSuffix = processedData.rxnSuffix;
            poolSuffix = processedData.poolSuffix;
            poolNames = unique(extractBefore(nameRxnsPool_i, rxnSuffix), 'stable');
            compNames = extractAfter(nameRxnsPool_i, rxnSuffix);
            compNames = unique(extractBetween(compNames, 1, 1), 'stable');
            compNames = [compNames; {'r'}];
            poolNames = strcat(poolNames, {'_'}, compNames);
            
            % separate indexes of rxns according to class:
            idxClasses = zeros(numel(idxRxnsPool_i),numel(poolNames));
            for j = 1:numel(poolNames)
                idxPool_j = contains(nameRxnsPool_i, poolNames{j});
                idxClasses(1:sum(idxPool_j),j) = idxRxnsPool_i(idxPool_j);
                poolNames{j} = [poolNames{j}, poolSuffix];
            end
            
        case {'AT1G73600'}
            idxClasses = idxRxnsPool_i(1);
            poolNames = cellstr(string(lipidPool_i));
            
            
        otherwise
            idxClasses = idxRxnsPool_i;
            poolNames = cellstr(string(lipidPool_i));
    end
    
    
    % (iv)add new pseudo-reaction:
    for k = 1:numel(poolNames)
        poolName_k = poolNames{k};
        idxRxnsPool_k = idxClasses(:,k);
        idxRxnsPool_k = idxRxnsPool_k(idxRxnsPool_k~=0);
    
        %lipidPool_id = processedData.poolName;
        model.S(size(model.S,1),end+1) = 0;
        model.rxns(end+1)           = {['pseudoRxn_',poolName_k]};
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

        % (v)add new pseudo-metabolite:
        model.S(end+1,idxRxnsPool_k) = 1;
        model.S(end,end)             = -1;
        model.b(end+1)               = 0;
        model.mets(end+1)            = {['pseudoMet_',poolName_k]};
        model.metNames(end+1)        = {'pseudo-metabolite'};
        model.csense(end+1)          = 'E';
        model.metCharges(end+1)      = 0;
        model.metFormulas(end+1)     = {''};
        model.metInChIString(end+1)  = {''};
        model.metKEGGID(end+1)       = {''};
        model.metChEBIID(end+1)      = {''};
        model.metPubChemID(end+1)    = {''};
        model.metLIPIDMAPSID(end+1)  = {''};


        % (vi) constrain the LB and UB of the pseudo-reaction:
        if nargin < 5
            flagUseRange = 0;
        end

        switch flagUseRange
            case 0
                ubPool = fsr_i*(fcPool+sdPool);

                if nargin < 6 % compute lb using fsr and relaxing with sd
                    lbPool = fsr_i*(fcPool-sdPool);

                else % relax the lb using the minimum value of the feasible space computed via FVA
                    minFluxPool = solFVA.minFlux(idxRxnsPool_k); % get minFlux for reactions
                    lbPool = min(minFluxPool);
                end

            case 1
                samplesWT = processedData.samplesWT;
                samples4FSR = samplesWT(idxRxnsPool_k,:);
                samples4FSR = sum(samples4FSR,1);
                meanSamples4FSR = mean(samples4FSR);
                lbPool = meanSamples4FSR*(fcPool-sdPool);
                ubPool = meanSamples4FSR*(fcPool+sdPool);

        end

        model.ub(end+1) = ubPool;
        model.lb(end+1) = lbPool;
    end

end

end

