function [ChainMets, BackboneMets] = obtainPseudometabolites(model)

% Chains list:
ChainMets = model.mets(contains(model.mets, '-chain') == 1);
ChainMets(contains(ChainMets,'Lipid-chain') ==1) = ''; % Eliminate Lipid-chain pseudometabolite from list
ChainMets = extractBefore(ChainMets, '[');

% Backbones list:
CompAbb = contains(model.compNames, {'cytosol','Cytosol','cytoplasm','Cytosol_0','cyt'}); % Find harmonized abbreviation for cytosol compartment
BackboneMets = model.mets(contains(model.mets, ['-backbone','[',model.comps{CompAbb},']']));
BackboneMets(contains(BackboneMets,'Lipid-backbone') ==1) = ''; % Eliminate Lipid-backbone pseudometabolite from lis
BackboneMets = extractBefore(BackboneMets, '[');

end