function model_mti = addConstrain2Pool(mt_i,mutantData,processedData_lusk,referenceFluxes,flagPools)
%##########################################################################
%
% This function is used to add extra constraints to pools of selected lipid
% species.
%
%##########################################################################


% get model previosuly created for mutant:
model_mti = processedData_lusk.modelsMT_w.(mt_i);
mutantData.model = model_mti;
mutantData.mt_i = mt_i;

switch flagPools
    case {'P'} % add constraints to pool of lipid rxn products
        % get reference fluxes for reaction products:
        referenceFluxesP = referenceFluxes.referenceFluxesP;
        
        % get standard deviation for each pool of metabolites:
        referenceFluxesP = getPoolSD(mt_i, referenceFluxesP, processedData_lusk);
        
        % add constraint to selected lipid pools:
        model_mti = computefluxsumPoolMT_mod(mutantData, referenceFluxesP, 0, [], 1);
        
    case {'R'} % add constraints to pool of lipid rxn reactants
        % get reference fluxes for reaction products:
        referenceFluxesR = referenceFluxes.referenceFluxesR;
        
        % get standard deviation for each pool of metabolites:
        referenceFluxesR = getPoolSD(mt_i, referenceFluxesR, processedData_lusk);
        
        % add constraint to selected lipid pools:
        model_mti = computefluxsumPoolMT_mod(mutantData, referenceFluxesR, 0, [], 1);

        
    case {'B'} % add constraints to pools of lipid rxn products and reactants
        referenceFluxesP = referenceFluxes.referenceFluxesP;
        referenceFluxesP = getPoolSD(mt_i, referenceFluxesP, processedData_lusk);
        model_mti = computefluxsumPoolMT_mod(mutantData, referenceFluxesP, 0, [], 1);

        mutantData.model = model_mti;
        referenceFluxesR = referenceFluxes.referenceFluxesR;
        referenceFluxesR = getPoolSD(mt_i, referenceFluxesR, processedData_lusk);
        model_mti = computefluxsumPoolMT_mod(mutantData, referenceFluxesR, 0, [], 1);
end

