function locusIDs = addLocusID(labelIDs, KOgenes, flagAverage)

%##########################################################################
%
% flagAverage   - (0) Default. Get locus ID for average lipid profiles
%                 (1) Get locus ID for each replicate of the lipid profiles
%
%##########################################################################

if nargin < 3
    flagAverage = 0;
end

switch flagAverage
    case 0
        locusIDs = repmat({'NA'}, 1, numel(labelIDs));
        count = 0;
        % get locus IDs:
        for kk = 1:numel(labelIDs)
            label_i = extractAfter(labelIDs{kk}, 'KO-');
            findID = strcmp(KOgenes(:,1), label_i);
            count = count + 1;

            if sum(findID)>0
                locusIDs{count} = KOgenes(findID,2);
            end
        end
        
    case 1
        locusIDs = repmat({'NA'}, 2, size(labelIDs,2));
        count = 0;
        
        % get information about replicates:
        replicateSuffix = extractAfter(labelIDs(2,:), '_');
        
        % get locus IDs:
        for kk = 1:size(labelIDs,2)
            label_i = extractAfter(labelIDs{1,kk}, 'KO-');
            findID = strcmp(KOgenes(:,1), label_i);
            count = count + 1;

            if sum(findID)>0
                locusIDs{2,count} = KOgenes(findID,2);
                
                % assign names to replicates:
                suffix_i = [KOgenes{findID,2}, '_' replicateSuffix{kk}];
                locusIDs{1,count} = suffix_i;
            end
        end
        
        % add names to Col-0 samples:
        idxCol0 = contains(labelIDs(1,:), 'Col0');
        locusIDs(2,idxCol0) = labelIDs(1, idxCol0);
        
        % add replicates info to Col-0 samples:
        locusIDs(1,idxCol0) = labelIDs(2, idxCol0);
        
end

locusIDs = cellstr(string(locusIDs));

