function [findLipid,matchedLipids] = CompareTailStructure(tailStructure, matchedMets, findMet, backCode)

% {'PS', 'PE', 'PC', 'PI', 'SQDG', 'PG','DG', 'MGDG', 'DGDG', 'TG', 'CL'} 
% Retrieve structural information from lipid:
tails = strsplit(tailStructure,{'/','_',';'});
tails{1,1} = strrep(tails{1,1}, '(', '');
tails{1,size(tails,2)} = strrep(tails{1,size(tails,2)}, ')', '');
tailSize = size(tails,2);

% Identify type of structural information available:
exactComposition = 0;
exactLocation = contains(tailStructure, '/');
randomLocation = contains(tailStructure, '_');
partialStructure = contains(tailStructure, ';');


% The data for the exact acyl chain composition is provided and the
% position of acyl chains is known
if sum(exactLocation) > 0
    exactComposition = 1;
    
% The data for the exact acyl chain composition is provided, but the
% position of acyl chains is unknown
elseif sum(randomLocation) > 0
    exactComposition = 2;

% The data for some of the acyl chains is available    
elseif sum(partialStructure) > 0
    exactComposition = 3;
end

if exactComposition == 1
    acylNumber = tailSize;
    [findMet, matchStructure] = CheckStereoLocation(matchedMets, tails, acylNumber, backCode, findMet);
    
elseif exactComposition == 3 || exactComposition == 0
    matchStructure = MatchTailsToMets(tails, matchedMets);
    
elseif exactComposition == 2
    switch tailSize
        case 1 % Data for just one (1) acyl chain is provided independently of the number of acyl chains in the molecule
            matchStructure = contains(matchedMets(:,2),tails(1,1));

        case 2 %Lipid species with 2 acyl chains
            matchStructure = contains(matchedMets(:,2),tails(1,1)).*contains(matchedMets(:,2),tails(1,2));

        case 3 %Lipid species with 3 acyl chains
            matchStructure = contains(matchedMets(:,2),tails(1,1)).*contains(matchedMets(:,2),tails(1,2))...
                .*contains(matchedMets(:,2),tails(1,3));

        case 4 %Lipid species with 4 acyl chains
            matchStructure = contains(matchedMets(:,2),tails(1,1)).*contains(matchedMets(:,2),tails(1,2))...
                .*contains(matchedMets(:,2),tails(1,3).*contains(matchedMets,tails(1,4)));
    end
end

findLipid = findMet(matchStructure > 0);
matchedLipids = matchedMets(matchStructure > 0,:);
        
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [findMet, matchStructure] = CheckStereoLocation(matchedMets, tails, acylNumber, backCode, findMet)

findTails = false(size(matchedMets,1), size(tails,2));
tailsModel = ObtainModelTails(backCode, matchedMets);

% First verify the metabolites in model with fixed stereospecific location:
for i = 1:size(tails,2)
    tailID = tails{i};
    getTailLocation = contains(tailsModel(:, i), tailID);
    findTails(:,i) = getTailLocation;
end

findTails = sum(findTails,2);
matchStructure = findTails == acylNumber;

% If none matching metabolites are found, check if metabolites in model
% have mixed stereospecific location:
if sum(matchStructure) == 0
    mixedStereoLocation = contains(matchedMets(:,2), '_');
    matchedMets = matchedMets(mixedStereoLocation, :);
    tailsModel = tailsModel(mixedStereoLocation, :);
    findMet = findMet(mixedStereoLocation);
    
    findTails = false(size(matchedMets,1), size(tails,2));
    for i = 1:size(tails,2)
        tailID = tails{i};
        for j = 1:size(tailsModel,1)
            getTailLocation = contains(tailsModel(j, :), tailID);
            findTails(j,i) = sum(getTailLocation,2);
        end
    end
    
    findTails = sum(findTails,2);
    matchStructure = findTails == acylNumber;
end

end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function matchStructure = MatchTailsToMets(tails, matchedMets)

matchStructure = false(size(matchedMets,1), size(tails,2));

for i = 1:size(tails,2)
    tailID = tails(1,i);
    findTail = contains(matchedMets(:,2),tailID);
    matchStructure(:,i) = findTail;
end

matchStructure = sum(matchStructure,2);

end