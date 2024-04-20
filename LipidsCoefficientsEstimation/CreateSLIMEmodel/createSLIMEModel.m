function [rxnName, SLIMErModel, biomassRxnsCE, LIPIDcoefficients] = createSLIMErModel(accessionData,Flag1,Flag2,biomassRxnsCE,varType)
%%-------------------------------------------------------------------------
% 'CreateSLIMErModel' function creates as an output a model that contains
% all the reactions of the Lipid Module and the corresponding
% stoichiometric coefficients of the SLIME pseudo-reactions, calculated
% using the experimental information provided. If required, the
% carboxylation/oxygenation ratio is also included.
% Hereinafter 'DD' refers to the number of days that the plants were
% subjected to the continuous darkness treatment.
%
% USAGE:
%   [rxnName, SLIMErModel, BiomassRxnsCE, LIPIDcoefficients] = createSLIMErModel(accessionData,Flag1,Flag2,BiomassRxnsCE)
%
% INPUTS:   
%   accessionData - structure containing the following fields:
%                   Day       - 0 (Control sample)
%                               3 (3DD)
%                               6 (6DD)
%                               6h  (Col-0, 6h)
%                               21h (Col-0, 21h)
%	
%                   Accession -	1           (when Day = 0)
%                               1 up to 282 (when Day = 3)
%                               1 up to 284 (when Day = 6)
%                               Col-0       (when Day = 6h or 21h)
%
%   BiomassRxnsCE - Structure containing the stoichiometric coefficients 
%                   for all the biomass components. This data has to be
%                   provided only in the case that biomass reactions are
%                   added to an already created model.
%
%       
%   Flag1         - 1 (Include carboxylation/oxygenation (C/O) ratio)
%                   0 (Do not include the C/O-ratio)
%
%   Flag2         - 1 (Create biomass reaction from scratch)
%                   0 (Add a new biomass reaction for accession_i under darkness)
%                   2 (Add a new biomass reaction for accession_i under light)
%
%   varType       - 1 (Data is retrieved for 1 accession)
%                 - 0 (Data is retrieved for multiple accessions)
%
% OUTPUTS:
%	SLIMErModel       - COBRA model structure that contains the merged list
%                       of metabolites from the Lipid module and the
%                       ModelTemplate, besides the stoichiometric 
%                       coefficients of the SLIME pseudo-reactions, 
%                       calculated using the experimental information 
%                       provided, and the accession- and condition-specific 
%                       biomass reaction, and the C/O-ratio. 
%
%   LIPIDcoefficients - Structure containing the stoichiometric
%                       coefficients for the lipid SLIME pseudo-reactions.
%
%   rxnName           - ID of the newly added accession- and condition-specific
%                       biomass reaction.
%--------------------------------------------------------------------------                       
%% 
% Initialize variables:
%           SpeciesName                SLIMEr_ID
poolNames  = {'CDP-DG'                'SLIMEr_146'
              'CL'                    'SLIMEr_328'
              'DG'                    'SLIMEr_411'
              'DGDG'                  'SLIMEr_456'
              'MGDG'                  'SLIMEr_798'
              'PA'                    'SLIMEr_870'
              'PC'                    'SLIMEr_908'
              'PE'                    'SLIMEr_1001'
              'PG'                    'SLIMEr_1068'
              'PGP'                   'SLIMEr_1117'
              'PI'                    'SLIMEr_1178'
              'PS'                    'SLIMEr_1271'
              'SQDG'                  'SLIMEr_1278'   
              'TG'                    'SLIMEr_1591'
              'Wax'                   'SLIMEr_1626'
              'Lipid-backbone[c]'     'SLIMEr_1627'
              'Lipid-chain[c]'        'SLIMEr_1628'}; 
      
if ~exist('Flag1', 'var')
    Flag1 = 1;
end

if ~exist('Flag2', 'var')
    Flag2 = 1;
end

if ~exist('varType', 'var')
    varType = 0;
end

if sum(varType) > 0
    accessionData.Accession = 1;
end

accession = accessionData.Accession;
day       = accessionData.Day;

switch Flag2
    case 1
        % Load model and obtain condition- and accesion-specific data:
        [biomassRxnsCE, model, LIPIDcoefficients] = createBiomassRxnsCE(accession,day);

        % Find index for carboxylation and oxygenation reactions:
        carbIdx = find(strcmp(model.rxns, 'RBC_h')); 
        oxyIdx  = find(strcmp(model.rxns, 'RBO_h'));

        % Create new metabolite: 'Biomass[c]'
        model = addMetabolite(model, 'Biomass[c]', 'metName','Biomass', 'Charge',0);

        % Add exchange reactions for missing Biomass metabolites:
        metsList    = {'dCTP[c]','dGTP[c]','dTTP[c]','Biomass[c]', 'cellulose2[c]'};
        LBs         = [0 0 0 0 0];
        UBs         = [1000 1000 1000 1000 1000];
        model       = addExchangeRxn(model, metsList, LBs, UBs);

        % Add Biomass[c] to metabolites list:
        biomassRxnsCE.CoeffsList = biomassRxnsCE.CoeffsList(:,:).*-1;
        biomassRxnsCE.MetsList   = [biomassRxnsCE.MetsList, {'Biomass[c]'}];
        coeffBiomass             = ones(1,size(biomassRxnsCE.CoeffsList,2));
        biomassRxnsCE.CoeffsList = [biomassRxnsCE.CoeffsList; coeffBiomass];

        % Adjust size of vectors:
        biomassRxnsCE.ChainsRxn     = [biomassRxnsCE.ChainsRxn; zeros(1,size(biomassRxnsCE.Accessions,2))];
        biomassRxnsCE.BackbonesRxn  = [biomassRxnsCE.BackbonesRxn; zeros(1,size(biomassRxnsCE.Accessions,2))];
        
        lipidsData = adjustSLIMErvector(poolNames,biomassRxnsCE);
    
    case 0
        poolNames  = [poolNames(:,1), strcat(poolNames(:,2),'_dark')];
        poolNames  = [strrep(poolNames(:,1),'[c]','[c2]'), poolNames(:,2)]; 
        lipidsData = adjustSLIMErvector(poolNames,biomassRxnsCE);
        model      = biomassRxnsCE.model;
        LIPIDcoefficients = '';

    case 2
        lipidsData = adjustSLIMErvector(poolNames,biomassRxnsCE);
        model      = biomassRxnsCE.model;
        LIPIDcoefficients = '';
end

% Add accession information:
switch varType
    case 1
        accessionData.accessionID = biomassRxnsCE.Accessions{varType};

    case 0
        accessionData.accessionID = biomassRxnsCE.Accessions{accession};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% STEP 1 %%%%% 
% Add biomass reaction for selected accession(s):
[rxnName, SLIMErModel] = addBiomassRxn(accessionData, biomassRxnsCE, model, lipidsData);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% STEP 2 %%%%% 
%%%%% Adjust ratio for RuBisCO carboxylation/oxygenation %%%%% 

% Set the RuBisCO ratio C/O according to measurements taken from 
% Tong et al. (2020) (https://doi.org/10.1038/s41467-020-16279-5)
if Flag1 == 1
    CO_ratio = 10.7/3.72;
    SLIMErModel.S(end+1,[carbIdx oxyIdx]) = [1 -CO_ratio];
    SLIMErModel.b(end+1) = 0;
    SLIMErModel.mets(end+1) = {'carb/oxy ratio'};
    SLIMErModel.metNames(end+1) = {'carb/oxy ratio'};

    % Complete information for remaining metabolites fields:
    SLIMErModel.csense          = [SLIMErModel.csense; 'E'];
    SLIMErModel.metCharges      = [SLIMErModel.metCharges; 0];
    SLIMErModel.metFormulas     = [SLIMErModel.metFormulas; {''}];
    SLIMErModel.metInChIString  = [SLIMErModel.metInChIString; {''}];
    SLIMErModel.metKEGGID       = [SLIMErModel.metKEGGID; {''}];
    SLIMErModel.metChEBIID      = [SLIMErModel.metChEBIID; {''}];
    SLIMErModel.metPubChemID    = [SLIMErModel.metPubChemID; {''}];
    SLIMErModel.metLIPIDMAPSID  = [SLIMErModel.metLIPIDMAPSID; {''}];
end

end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function lipidsData = adjustSLIMErvector(poolNames,biomassRxnsCE)

% Adjust SLIMEr coefficients vector in model.S:
idxSubspecie            = ismember(poolNames(:,1), biomassRxnsCE.subSpeciesRxn.Names);
lipidsData.subspecieIdx = poolNames(idxSubspecie,:);
lipidsData.pseudoLipid  = poolNames(16:17,:);

end