function newDataSet = removeIsNaN(dataSet, flag)

if ~exist('flag','var')
    flag = 0;
end

dataSet(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),dataSet)) = {''};
idx = cellfun(@ischar, dataSet);
dataSet(idx) = cellfun(@(x) string(x), dataSet(idx), 'UniformOutput', false);
dataSet = cellstr(string(dataSet));

% Replace non-importable data:
idx = strcmpi(string(dataSet),"ActiveX VT_ERROR: ");
dataSet(idx) = {''};

newDataSet = dataSet;

switch flag
    case 1
        newDataSet = newDataSet(~cellfun(@isempty, newDataSet(:,2)),:);
        
        % Exclude metabolites different to lipids
        newDataSet = [newDataSet(:,1),newDataSet(:,54:end)];
        
        % Adjust names of metabolites:
        newDataSet = adjustNames(newDataSet);

    case 2
        % Adjust names of metabolites:
        newDataSet = adjustNames(newDataSet);      
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Adjust name of lipid species to eliminate reference to isomers:
function newDataSet = adjustNames(newDataSet)

for i = 1:size(newDataSet,2)
    lipidName = newDataSet{1,i};
    if contains(lipidName, '(')
        lipidName = extractBefore(newDataSet{1,i}, ' (');
        newDataSet{1,i} = lipidName;
    end
end

end
