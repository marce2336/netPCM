function coeffsVector = changeRxnCoefficient(data, model, coeffsVector, accession, flag)
%%-------------------------------------------------------------------------
% 'ChangeRxnCoefficient' function adjusts the coefficients vector for the
% metabolites paricipating in the SLIME pseudo-reactions to match the
% dimensions of model.S.
%
% USAGE:
%   coeffsVector = changeRxnCoefficient(data, model, coeffsVector, accession, flag)
%
% INPUTS:   
%
%	data            -   Structure containing the IDs and abundances of
%                       metabolites paricipating in the SLIME
%                       pseudo-reactions.
%
%   model           -   COBRA model structure.
%
%   coeffsVector    -   Coefficients vector that is going to be adjusted.
%
%   accession       -   ID of accession for which the coefficients vector
%                       is adjusted.
%
%   flag            -   Environmental conditions for which the experimental
%                       data is provided:
%                       [] - Control conditions
%                       1  - Light conditions
%                       2  - Dark conditions
%
%
% OUTPUTS:
%	coeffsVector   -    Adjusted coefficients vector for the metabolites 
%                       paricipating in the SLIME pseudo-reactions.
%--------------------------------------------------------------------------
%%
if ~exist('flag', 'var')
    abbreviation = '[c]';
    
elseif flag == 1 % 1 stands for light conditions
    abbreviation = '[c1]';
    
elseif flag == 2 % 2 stands for dark conditions
    abbreviation = '[c2]';
end

for j = 1:size(data.dataNames,1) %Chain
    getMetID = [data.dataNames{j,1},abbreviation];
    metIdx = strcmp(model.mets, getMetID);
    coeffsVector(metIdx,accession) = data.dataAbundance(j,accession);
end

end