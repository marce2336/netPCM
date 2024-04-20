function findIdxs = getSubspeciesIdx(poolName_i, backbones)
%%-------------------------------------------------------------------------
% 'GetSubspeciesIdx' function finds the indexes for lipid subespecies
% within each lipid species pool.
%
% USAGE:
%   findIdxs = GetSubspeciesIdx(poolName_i, backbones)
%
% INPUTS:   
%   poolName_i  -   ID of lipid pool.
%
%   backbones	-	List with IDs for all lipid backbones of the Lipid
%   Module.
%
% OUTPUTS:
%	findIdxs	-	List with indexes for lipid subespecies for a
%                   determined lipid species pool.
%                                  
%--------------------------------------------------------------------------
%%
pattern    = [poolName_i '+[0-9]' '+[0-9]'];
matchSubspecie = regexp(backbones, pattern, 'match');
findIdxs = (~cellfun(@isempty, matchSubspecie)).*(strlength(poolName_i)+2 == strlength(backbones)).*(contains(backbones, poolName_i)) == 1;
findIdxs = find(findIdxs);

end
