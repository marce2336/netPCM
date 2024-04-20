function [pFBA_sol, bio_opt] = getFluxDist_pFBA(model, solver)
%%-------------------------------------------------------------------------
% 'GetFluxDist_pFBA' function was created for estimating the reference flux
% distributions via pFBA for Arabidopsis Col-0.
% This code was modified from the original published by Tong, Küken &
% Nikoloski (2020), https://doi.org/10.1038/s41467-020-16279-5.
%
% USAGE:
%   [pFBA_sol, bio_opt] = GetFluxDist_pFBA(model)
%
% INPUTS:   
%   model         - COBRA model structure.
%
% OUTPUTS:
%	pFBA_sol      - Structure containing the solution object:
%                   f - Objective value
%                   v - Reaction rates
%                   y - Dual for the metabolites
%                   w - Reduced costs of the reactions
%                   s - Slacks of the metabolites
%                   stat - Solver status in standardized form:
%                          -1 - No solution reported (timelimit, numerical problem etc)
%                           1 - Optimal solution
%                           2 - Unbounded solution
%                           0 - Infeasible                
%
%   bio_opt       - ID of the newly added accession- and condition-specific
%                   biomass reaction.
%--------------------------------------------------------------------------
%%
% Set 'gurobi' as default solver if not specified:
if ~exist('solver', 'var')
    solver = 'gurobi';
end

%% Calculate Arabidopsis Col-0 flux distribution via pFBA %%%%% 
    
    % Run FBA after adjusting C/O ratio:
    changeCobraSolver(solver, 'LP', 0);
    sol = optimizeCbModel(model);
    bio_opt = sol.f; % biomass under optimal growth condition and C/O ratio
    
%% Run pFBA
    % At optimum biomass minimize flux through gene associated reactions
    model.lb(model.c~=0) = bio_opt;
    model.ub(model.c~=0) = bio_opt;

    gene_associated_rxns = full(sum(model.rxnGeneMat,2)) > 0;
    
    %gene_associated_rxns=find(sum(model.rxnGeneMat(:,find(cellfun(@isempty,strfind(model.genes,'AT'))==0))')>0);
    model.c=zeros(size(model.c));
    model.c(gene_associated_rxns)= -1;

    pFBA_sol = optimizeCbModel(model);
    
end

