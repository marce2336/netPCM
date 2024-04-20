function [biomassRxnsCE, model, LIPIDcoefficients] = createBiomassRxnsCE(accession,day)
%%-------------------------------------------------------------------------
% 'CreateBiomassRxnsCE' function generates the stoichiometric coefficients
% for all the biomass components for the selected accession and day, using
% the experimental data provided.
% Hereinafter 'DD' refers to the number of days that the plants were
% subjected to the continuous darkness treatment.
%
% USAGE:
%   [biomassRxnsCE, model, LIPIDcoefficients] = createBiomassRxnsCE(accession,day)
%
% INPUTS:   
%   day       - 0   (Control sample)
%               3   (3DD)
%               6   (6DD)
%               6h  (Col-0, 6h)
%               21h (Col-0, 21h)
%	
%   accession -	1           (when Day = 0)
%               1 up to 282 (when Day = 3)
%               1 up to 284 (when Day = 6)
%               Col-0       (when Day = 6h or 21h)
%               All         (when it is necessary to calculate the
%                               stoichiometric coefficients for all the 
%                               accessions in 3DD and 6DD)
%       
% OUTPUTS:
%	model             - COBRA model structure that contains the merged list
%                       of metabolites from the Lipid module and the
%                       ModelTemplate.
%
%   biomassRxnsCE     - Structure containing the stoichiometric
%                       coefficients for all the biomass components.
%
%   LIPIDcoefficients - Structure containing the stoichiometric
%                       coefficients for the lipid SLIME pseudo-reactions.
%--------------------------------------------------------------------------
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1). Obtain coefficients for biomass components for selected accession and
% treatment:
[biomassCoeff, model, LIPIDcoefficients] = getBiomassComponents(accession,day);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2). Extract data for lipid species (backbones) and subspecies:
%Initialize variables:
poolNames  = {'CDP-DG', 'CL', 'DG', 'DGDG', 'MGDG', 'PA', 'PC', 'PE',...
    'PG', 'PGP', 'PI', 'PS', 'SQDG', 'TG', 'Wax'};
backbonesSpecies   = LIPIDcoefficients.BackboneNames;
backbonesSpecies   = extractBefore(backbonesSpecies, '-backbone');
subIdxs            = zeros(size(LIPIDcoefficients.BackboneNames,1),1);
count              = 1;

% Find  matching subspecies for each pool:
for i = 1:size(poolNames,2)
    poolName_i = poolNames{i};
    findIdxs   = getSubspeciesIdx(poolName_i, backbonesSpecies);
    subIdxs(count:count+(size(findIdxs,1))-1) = findIdxs;
    count      = count + size(findIdxs,1);
end

% Create list with subspecies names and coefficients:
subIdxs         = subIdxs(subIdxs ~= 0);
subSpeciesNames = LIPIDcoefficients.BackboneNames(subIdxs);
LIPIDcoefficients.subSpeciesNames = subSpeciesNames;
subSpeciesCoeff = LIPIDcoefficients.BackboneAbundance(subIdxs,:);
LIPIDcoefficients.subSpeciesCoeff = subSpeciesCoeff;
    
% Eliminate subspecies from species (backbones) list:
LIPIDcoefficients.BackboneNames(subIdxs) = '';
LIPIDcoefficients.BackboneAbundance(subIdxs,:) = '';  
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3). Re-scale coefficients of chains, backbones and subspecies:
for i = 1:size(LIPIDcoefficients.Accessions,2)
    % Re-scale chains:
    totalChains  = sum(LIPIDcoefficients.ChainAbundance(:,i));
    scaledChains = LIPIDcoefficients.ChainAbundance(:,i)./totalChains;
    LIPIDcoefficients.ChainAbundance(:,i) = scaledChains;
    
    % Re-scale subspecies belonging to the same species (backbone) pool:
    backbonesSE   = LIPIDcoefficients.subSpeciesNames;
    backbonesSE   = extractBefore(backbonesSE, '-backbone');
    for j = 1:size(poolNames,2)
        poolName_i   = poolNames{j};
        findIdxs     = getSubspeciesIdx(poolName_i, backbonesSE);
        totalPool    = sum(LIPIDcoefficients.subSpeciesCoeff(findIdxs,i));
        
        % Add total pool abundance to backbones list:
        idxBackbone = startsWith(LIPIDcoefficients.BackboneNames, [poolName_i,'-backbone']);
        LIPIDcoefficients.BackboneAbundance(idxBackbone,i) = totalPool;
        
        % Rescale coefficients for subspecies:
        reScaledPool = LIPIDcoefficients.subSpeciesCoeff(findIdxs,i)./totalPool;
        LIPIDcoefficients.subSpeciesCoeff(findIdxs,i) = reScaledPool;
    end
    
    % Re-scale backbones:
    totalBackbones  = sum(LIPIDcoefficients.BackboneAbundance(:,i));
    scaledBackbones = LIPIDcoefficients.BackboneAbundance(:,i)./totalBackbones;
    LIPIDcoefficients.BackboneAbundance(:,i) = scaledBackbones;
    
end

LIPIDcoefficients.subSpeciesCoeff(isnan(LIPIDcoefficients.subSpeciesCoeff)) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4). Create coefficient vectors for lipid chains in accession-specific SLIMEr:

% Initialize variables:
biomassRxnsCE.Accessions = LIPIDcoefficients.Accessions;
biomassRxnsCE.MetsList   = biomassCoeff.Mets';
biomassRxnsCE.CoeffsList = biomassCoeff.S;
biomassRxnsCE.LBs        = -1000.*ones(length(biomassCoeff.Mets),1);
biomassRxnsCE.UBs        = 1000.*ones(length(biomassCoeff.Mets),1);

% Adjust sign of coefficients for backbones and chains being consumed in
% the reaction:
LIPIDcoefficients.BackboneAbundance = -1.*LIPIDcoefficients.BackboneAbundance;
LIPIDcoefficients.ChainAbundance    = -1.*LIPIDcoefficients.ChainAbundance;
LIPIDcoefficients.subSpeciesCoeff   = -1.*LIPIDcoefficients.subSpeciesCoeff;

% Add Lipid-backbone and Lipid-chain pseudometabolites to list:
LIPIDcoefficients.BackboneNames     = [LIPIDcoefficients.BackboneNames; 'Lipid-backbone'];
LIPIDcoefficients.BackboneAbundance = [LIPIDcoefficients.BackboneAbundance;...
    ones(1,size(LIPIDcoefficients.BackboneAbundance,2))];

LIPIDcoefficients.ChainNames     = [LIPIDcoefficients.ChainNames; 'Lipid-chain'];
LIPIDcoefficients.ChainAbundance = [LIPIDcoefficients.ChainAbundance;...
    ones(1,size(LIPIDcoefficients.ChainAbundance,2))];

% Adjust coefficients of the reaction vector for Lipid chains-pseudoreactions:
chainsRxn            = zeros(size(model.S,1),size(LIPIDcoefficients.Accessions,2));
chains.dataNames     = LIPIDcoefficients.ChainNames;
chains.dataAbundance = LIPIDcoefficients.ChainAbundance;

for i = 1:size(LIPIDcoefficients.Accessions,2)
    chainsRxn = changeRxnCoefficient(chains, model, chainsRxn, i);
end
biomassRxnsCE.ChainsRxn = chainsRxn;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5). Create coefficient vectors for lipid subspecies in accession-specific SLIMEr:

% Create an accession specific stoichiometric coefficients vector for each subespecies:
subSpeciesRxn.Names = cell(size(poolNames,2),1);
count = 0;
for i = 1:size(poolNames,2)
    % First find out if the pool has a non-zero coefficient:
    poolName_i = poolNames{i};
    findPool_i = strcmp(LIPIDcoefficients.BackboneNames, [poolName_i '-backbone']);
    coeffPool_i = LIPIDcoefficients.BackboneAbundance(findPool_i);
    
    % If species has a non-zero value proceed to modify the coefficients for subspecies:
    if sum(coeffPool_i) ~= 0
        count = count + 1;
        subSpeciesRxn.Names{count} = poolName_i;
        
        % Find subespecies from pool:
        findIdxs = getSubspeciesIdx(poolName_i, backbonesSE);
        subSpecies.dataNames = [LIPIDcoefficients.subSpeciesNames(findIdxs); [poolName_i '-backbone']];
        subSpecies.dataAbundance = LIPIDcoefficients.subSpeciesCoeff(findIdxs,:);
        subSpecies.dataAbundance = [subSpecies.dataAbundance; ones(1, size(subSpecies.dataAbundance,2))];
        coeffsSubspecie_i = zeros(size(model.S,1),size(LIPIDcoefficients.Accessions,2));
        
        for j = 1:size(LIPIDcoefficients.Accessions,2)
            coeffsSubspecie_i = changeRxnCoefficient(subSpecies, model, coeffsSubspecie_i, j);
        end
        
        subSpeciesRxn.(poolName_i) = coeffsSubspecie_i;
    end
end

subSpeciesRxn.Names = subSpeciesRxn.Names(~cellfun(@isempty, subSpeciesRxn.Names));
biomassRxnsCE.subSpeciesRxn = subSpeciesRxn;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6). Create coefficient vectors for lipid backbones in accession-specific SLIMEr:

% Adjust coefficients of the reaction vector for Lipid backbones-pseudoreactions:
backbonesRxn            = zeros(size(model.S,1),size(LIPIDcoefficients.Accessions,2));
backbones.dataNames     = LIPIDcoefficients.BackboneNames;
backbones.dataAbundance = LIPIDcoefficients.BackboneAbundance;

for i = 1:size(LIPIDcoefficients.Accessions,2)
    backbonesRxn = changeRxnCoefficient(backbones, model, backbonesRxn, i);
end

biomassRxnsCE.BackbonesRxn = backbonesRxn;

end
