function [matchedMets, findMet] = MatchSumComposition(Hits, model, tailCode, backCode)

% Get metabolite sum composition:
tail = strrep(tailCode, ':', '_');

% Get names of matched metabolites:
matchedMets = model.metNames(Hits);

% Get abbreviations of matched metabolites:
matchedAbbs = model.mets(Hits);

% Keep only lipids with the same sum composition (mass isomers):
matchedSumComposition = contains(matchedAbbs, tail);
matchedMets = matchedMets(matchedSumComposition);

% Get indexes of mass isomers:
findMet = find(Hits);
findMet = findMet(matchedSumComposition);

% Split up lipids into their components:
if size(matchedMets,1) == 1
    matchedMets = cellstr(strsplit(string(matchedMets), backCode));  
elseif size(matchedMets,1) > 1
    matchedMets = split(matchedMets, backCode);
end

end