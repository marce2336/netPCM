function [rxnName, SLIMErModel] = addBiomassRxn(accessionData, biomassRxnsCE, model, lipidsData)
%%-------------------------------------------------------------------------
% 'addBiomassRxn' function creates as an output a model that contains
% all the reactions of the Lipid Module and the corresponding
% stoichiometric coefficients of the SLIME pseudo-reactions, calculated
% using the experimental information provided.
% Hereinafter 'DD' refers to the number of days that the plants were
% subjected to the continuous darkness treatment.
%
% USAGE:
%   [rxnName, SLIMErModel] = addBiomassRxn(accessionData, biomassRxnsCE, model, lipidsData)
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
%   BiomassRxnsCE - Structure containing the stoichiometric coefficients
%                   for all the biomass components.
%
%   model         - COBRA model structure that contains the merged list of
%                   metabolites from the Lipid module and the ModelTemplate
%
%   lipidsData    - Structure containing the information about the lipid
%                   pseudo-species used in the SLIME pseudo-reactions. 
%
% OUTPUTS:
%	SLIMErModel   - COBRA model structure that contains the merged list of
%                   metabolites from the Lipid module and the
%                   ModelTemplate, besides the stoichiometric coefficients
%                   of the SLIME pseudo-reactions, calculated using the
%                   experimental information provided, and the
%                   accession- and condition-specific biomass reaction.
%
%   rxnName       - ID of the newly added accession- and condition-specific
%                   biomass reaction.
%--------------------------------------------------------------------------
%%
subSpecie = lipidsData.subspecieIdx;
pseudoLipid = lipidsData.pseudoLipid;
accessionID = accessionData.accessionID;
accession_i = accessionData.Accession;

model = addSLIMErCoefficients(model, subSpecie, biomassRxnsCE, pseudoLipid, accession_i);

if isfield(biomassRxnsCE,'dark') && strcmp(biomassRxnsCE.dark, 'true')
    metsList = strrep(biomassRxnsCE.MetsList, ']', '2]');
else
    metsList = biomassRxnsCE.MetsList;
end

% Create Accession- and Condition-specific biomass reaction:
rxnName = ['Bio_',accessionID];
model   = addReaction(model, rxnName, 'metaboliteList', metsList,...
'stoichCoeffList', biomassRxnsCE.CoeffsList(:,accession_i), 'reversible', false);

% Change objective function:
model = changeObjective(model, rxnName, 1);

SLIMErModel = model;    

end
