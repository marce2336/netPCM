function model = addSLIMErCoefficients(model, isSubspecie, biomassRxnsCE, pseudoLipid, accession_i)
%%-------------------------------------------------------------------------
% 'addSLIMErCoefficients' function adds the stoichiometric coefficients to
% the SLIME pseudo-reactions using the experimental information provided,
% to a defined accession under the selected conditions.
%
% USAGE:
%   model = addSLIMErCoefficients(model, isSubspecie, BiomassRxnsCE, pseudoLipid, accession_i)
%
% INPUTS:   
%   accessionData - structure containing the following fields:
%                   Day       - 0 (Control sample)
%                               3 (3DD)
%                               6 (6DD)
%	
%                   Accession -	1           (when Day = 0)
%                               1 up to 282 (when Day = 3)
%                               1 up to 284 (when Day = 6)
%       
%   model         - COBRA model structure that contains the merged list of
%                   metabolites from the Lipid module and the ModelTemplate
%
%   isSubspecie   - List with lipid subspecies available in the 
%                   experimental data.
%
%   BiomassRxnsCE - Structure containing the stoichiometric coefficients
%                   for all the biomass components.
%
%   pseudoLipid   - List with IDs for metabolites participating in the
%                   SLIME pseudo-reactions.
%
%   accession_i   - ID of the accession for which the biomass reaction is
%                   being created.
%
%
% OUTPUTS:
%	model         - COBRA model structure that contains the merged list of
%                   metabolites from the Lipid module and the
%                   ModelTemplate, besides the stoichiometric coefficients
%                   of the SLIME pseudo-reactions, calculated using the
%                   experimental information provided.
%%--------------------------------------------------------------------------
%%

for i = 1:size(isSubspecie,1)
    subspecie_i = isSubspecie{i,1};
    SLIMEr_i = isSubspecie{i,2};
    findRxnID = strcmp(model.rxns, SLIMEr_i);
    
    if (isfield(biomassRxnsCE,'dark') && strcmp(biomassRxnsCE.dark, 'true')) || (isfield(biomassRxnsCE,'light') && strcmp(biomassRxnsCE.light, 'true'))
        model.S(:,findRxnID) = biomassRxnsCE.subSpeciesRxn.(subspecie_i)(:,accession_i);
        
    else
        model.S(:,(findRxnID==1)) = [biomassRxnsCE.subSpeciesRxn.(subspecie_i)(:,accession_i); 0];
    end
end
              
% Adjust Backbones-pseudoreaction coefficients
for i = 1:size(pseudoLipid, 1)
    subspecie_i = pseudoLipid{i,1};
    SLIMEr_i = pseudoLipid{i,2};
    getRxnID = strcmp(model.rxns, SLIMEr_i);
    
    lipidSLIMEr = 0;
    if sum(strcmp(subspecie_i, 'Lipid-backbone[c]')) > 0 || sum(strcmp(subspecie_i, 'Lipid-backbone[c2]')) > 0
        lipidSLIMEr = 1;
    end
    
    switch lipidSLIMEr
        case 1
            model.S(:,getRxnID) = biomassRxnsCE.BackbonesRxn(:,accession_i);
        case 0
            model.S(:,getRxnID) = biomassRxnsCE.ChainsRxn(:,accession_i);
    end
end

end