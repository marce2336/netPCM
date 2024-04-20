function [growthMTs, fba_KOs, modelsMTs] = optimizeGrowthMT(uniqueListKOs, processedData, flagOpt, flagSolver)
%##########################################################################
%
% This function optimizes growth of mutants (MTs) by considering three scenarios:
% (1) Optimization under constraints imposed on lipid-related reactions
%     using flux sum computed for mutants.
%
% (2) Optimization under constraints imposed on lipid-related reactions
%     using flux sum computed for mutants, besides setting to zero the
%     fluxes of reactions catalyzed by KO gene product.
%
% (3) Optimization under constraints imposed on lipid-related reactions
%     using flux sum computed for mutants, besides setting to zero the
%     fluxes of reactions catalyzed by KO gene product, but only for
%     reactions catalyzed by unique enzymes.
%
% flagOpt  -  Default. (1) Optimize growth under constraints of case scenario (1)
%             (2) Optimize growth under constraints of case scenario (2)
%             (3) Optimize growth under constraints of case scenario (3)
%
% flagRelaxation - (true) Use sd to calculate lower and upper bounds
%                  (false) Default. Do not relax upper and lower bounds
%
% flagSolver - ('gurobi') 
%              ('ibm_cplex') Default
%
%##########################################################################

if nargin < 3
    flagOpt = 1;
    flagSolver = 'ibm_cplex';
end


% Change solver:
if nargin < 4
    flagSolver = 'ibm_cplex';
end

changeCobraSolver(flagSolver, 'LP')


% Pre-allocate space to save predicted growth MTs:
growthMTs = cell(size(uniqueListKOs,1),2);

% create structs for saving the necessary information:
processedProfiles = processedData.processedProfiles;


% compute sum of incoming fluxes for each of the lipid pools:
referenceFluxes_wSFR = processedData.referenceFluxes;
referenceFluxes_wSFR = sumIncomingFluxes(referenceFluxes_wSFR);


% Optimize growth for MTs taking into account case scenario selected:
for i = 2:size(uniqueListKOs,1)
    ko_i = uniqueListKOs{i,2};
    growthMTs{i,1} = ko_i;

    % get lipid profile for mutant:
    idxKO_i = strcmp(processedProfiles.FCaverages(1,:), ko_i);

    if sum(idxKO_i)>0
        idxKO_i(1) = true;

        % constrain model using the flux sum calculated for the mutants:
        modelMT = createModelMT(processedData, idxKO_i, referenceFluxes_wSFR);

        % add additional constraints if this option is selected:
        switch flagOpt
            case 1
                growthMTs(1,:) = {'mutant', 'growth_case1'};

            case 2
                growthMTs(1,:) = {'mutant', 'growth_case2'};

                % find reaction(s) catalyzed by KO gene:
                findRxnKO = contains(modelMT.grRules, ko_i);

                % impose zero flux for the knocked out reaction(s):
                modelMT.ub(findRxnKO) = 0;

            case 3
                growthMTs(1,:) = {'mutant', 'growth_case3'};

                % find reaction(s) catalyzed by KO gene:
                findRxnKO = contains(modelMT.grRules, ko_i);

                % identify if any of the rxn(s) is catalyzed by isoenzyme:
                findIso = contains(modelMT.grRules, ko_i).*(contains(modelMT.grRules, 'or')) == 1;

                % impose zero flux for the knock-out reaction(s), only
                % for rxns catalyzed by unique gene product:
                findRxnKO(findIso) = false;
                modelMT.ub(findRxnKO) = 0;

        end
    
        % optimize growth of mutant:
        fba_sol = optimizeCbModel_modified(modelMT);
        growthMTs{i,2} = cellstr(string(fba_sol.f));
    
        fba_KOs.(ko_i).v = fba_sol.v;
        fba_KOs.(ko_i).f = fba_sol.f;
        
        % save MT model for further analyzes:
        modelsMTs.(ko_i) = modelMT;
    end
end


