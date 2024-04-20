%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data = GetSLIMErCoeff(data_old,model,condense)
%
% Modified from: Benjamin J. Sanchez (2018-11-18)
% data = convertEjsingData(data,model,condense)
% 
% The units of the input data must be provided in g g-1DW
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = GetSLIMErCoeff(data_old,model,condense)

%Map each lipid to corresponding positions in the model:
MWs        = zeros(size(data_old.metNames));
backbones  = cell(size(data_old.metNames));
ChainsAll  = cell(size(data_old.metNames));
snPositions = zeros(size(data_old.metNames));

for i = 1:length(data_old.metNames)
    [pos, backbone, chains] = matchToModel(model,data_old.metNames{i});
    if sum(pos) > 0
        MWi            = getMWfromFormula(model.metFormulas(pos));% MW units: g/mmol
        MWs(i)         = mean(MWi); % MW units: g/mmol
        backbones{i}   = backbone;
        ChainsAll{i}   = chains;
        snPositions(i) = size(cellstr(chains),2); % Number of acyl chains in lipid
    end
end

%Filter out mets with no MW computed (i.e. not in model):
backbones          = backbones(MWs > 0);
ChainsAll          = ChainsAll(MWs > 0); 
data_old.metNames  = data_old.metNames(MWs > 0);
data_old.abundance = data_old.abundance(MWs > 0); % Data is provided in g g-1DW
data_old.std       = data_old.std(MWs > 0); % Data is provided in g g-1DW
snPositions        = snPositions(MWs > 0);
MWs                = MWs(MWs > 0);

if condense
    %Initialize variables
    data.lipidData.metNames     = data_old.metNames;
    data.lipidData.abundance    = data_old.abundance;
    data.lipidData.std          = data_old.std;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Backbones: get backbones abundance
    %Change units to mmol g-1DW and then add up
    data.lipidData              = changeUnits(data.lipidData,MWs);
    abundanceData               = data.lipidData;
    data.lipidData              = squashLipids(data.lipidData,backbones);
    data.lipidData.Accession    = data_old.Accession;
    
    %Estimate total lipid content in g g-1DW and mmol g-1DW:
    data.lipidData.TotalLipid_g     = sum(data_old.abundance);
    data.lipidData.TotalLipid_umol  = sum(data.lipidData.abundance);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Tails: First get tail names listed in the model
    pathMets = fullfile('InputData','AcylChainsList.xlsx');
    [~, ~, chainNames] = xlsread(pathMets,'AcylChainsList','A2:C85');
    chainNames(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),chainNames)) = {''};
    idx = cellfun(@ischar, chainNames);
    chainNames(idx) = cellfun(@(x) string(x), chainNames(idx), 'UniformOutput', false);
    chainNames = cellstr(string(chainNames));
    clearvars pathMets idx
    
    %Harmonize tail names:
    tempNames     = strrep(chainNames,'-chain[c]','');
    tempNames     = strrep(tempNames,'_',':');
    
    for i = 1:size(tempNames,1)
        chainID = tempNames{i,1};
        if strlength(chainID) <= 6 && ~contains(chainID,'h')
            tempNames{i,1}     = strrep(chainID,'C','');
        end
    end
    
    %Now add up:
    data.chainData.metStructure = ChainsAll;
    data.chainData.abundance    = abundanceData.abundance; % Units: umol g-1DW
    data.chainData.std          = abundanceData.std;
    data.chainData.snPositions  = snPositions;
    data.chainData              = squashLipids(data.chainData,tempNames);
    
    %Correct fields by adding chain ID according to model nomenclature:
    for iii = 1:length(data.chainData.metNames)
        chainID = data.chainData.metNames{iii,1};
        if strlength(chainID) <= 4
            data.chainData.metNames{iii} = ['C',chainID,'-chain'];
        else
            data.chainData.metNames{iii} = [chainID,'-chain'];
        end
    end
    
    data.chainData.metNames = strrep(data.chainData.metNames, ':','_');
    
else
    %Only change units to g/gDW:
    %data     = changeUnits(data_old,MWs);
    %data.MWs = MWs;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function data = squashLipids(data_old,metNames)

%Initialize variables:
data.metNames  = unique(metNames(:,1));
data.abundance = zeros(size(data.metNames));
data.std       = zeros(size(data.metNames));

%Determine if backbones or tails are being squashed:
if length(data.metNames) == length(metNames)
    isTail = true;
    % Find number of mass isomers in the metabolite pool:
    getMassIsomers = cellfun('size',data_old.metStructure, 1);

    % Get total number of acyl chains in the metabolite pool:
    totalFAs = data_old.snPositions.*getMassIsomers;

    % Calculate abundance of each acyl chain in the metabolite pool:
    abundanceFAs = (data_old.abundance.*data_old.snPositions)./totalFAs;

else
    isTail = false;
end
        
%Squash each species by adding all abundances and averaging std:
for i = 1:length(data.metNames)
    % Squash chains:
    if isTail
         % Find chain in metabolite pools:
        getChain    = data.metNames(i);
        
        % Find abundances for mass isomers:
        isChain1 = cellfun(@(x) strcmp(x(:), getChain), data_old.metStructure, 'UniformOutput', false);
        hits1 = cellfun(@sum, isChain1);
        
        % Find abundances for pools with unique species:
        isChain2 = double(strcmp(data_old.metStructure, getChain));
        
        % Estimate total abundance of acyl chain:
        data.abundance(i) = sum(abundanceFAs.*hits1) + sum(abundanceFAs.*isChain2);
        data.std(i)       = mean((data_old.std.*hits1) + (data_old.std.*isChain2));
        
    % Squash backbones:    
    else
        hits = strcmp(metNames,data.metNames{i});
        data.abundance(i) = sum(data_old.abundance.*hits);
        data.std(i)       = mean(data_old.std.*hits);
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function data = changeUnits(data_old,MWs)

%Transform units from g g-1DW to mmol g-1DW:
data.metNames   = data_old.metNames;
data.abundance  = data_old.abundance./MWs.*1e+03;   %(mmol lipid i)/gDW
data.std        = data_old.std./MWs.*1e+03;         %(mmol(lipid i)/gDW

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%