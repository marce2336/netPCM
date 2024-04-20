function selectedKOs = findKOmodelGenes(flagSave, flagFilter)
%##########################################################################
%
% This function identifies the model genes for which T-DNA lines are available:
% flagSave - (1) save output variable
%            (0) don't save
%
% flagFilter - (1) Filter data to leave only mutants for genes included in model
%              (0) Use the complete dataset
%##########################################################################

if ~exist('flagSave', 'var')
    flagSave = 0;
end

if ~exist('flagFilter', 'var')
    flagFilter = 1;
end

% Load list of T-DNA lines for which the lipidomics data was measured and
% obtain the numbers used to identify the samples:
pathFile = fullfile('InputData', 'sample_description_RunningOrder.txt');

opts = delimitedTextImportOptions("NumVariables", 4);
opts.DataLines = [2, Inf];
opts.Delimiter = "\t";
opts.VariableNames = ["Samplename", "Measurementday", "Tray", "Runningorder"];
opts.VariableTypes = ["char", "double", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts = setvaropts(opts, "Samplename", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "Samplename", "EmptyFieldRule", "auto");
listKOsMeasured = readtable(pathFile, opts);
listKOsMeasured = table2cell(listKOsMeasured);
numIdx = cellfun(@(x) ~isnan(str2double(x)), listKOsMeasured);
listKOsMeasured(numIdx) = cellfun(@(x) {str2double(x)}, listKOsMeasured(numIdx));
listKOsMeasured = [opts.SelectedVariableNames; listKOsMeasured];
listKOsMeasured = cellstr(string(listKOsMeasured));

sampleNumbers = listKOsMeasured(2:end, 1);
sampleNumbers = extractBefore(sampleNumbers, '_');
sampleNumbers = unique(sampleNumbers);
sampleNumbers(contains(sampleNumbers, 'T')) = ''; % Eliminate control samples from list

% Load list of genes identity:
pathFile = fullfile('InputData', 'mutant_labels_and_descriptions.txt');
opts.VariableNames = ["Tube", "number Locus", "Gene", "Gene function"];
opts.VariableTypes = ["double", "char",  "char", "char"];
listKOGenes = readtable(pathFile, opts);
listKOGenes = table2cell(listKOGenes);
numIdx = cellfun(@(x) ~isnan(str2double(x)), listKOGenes);
listKOGenes(numIdx) = cellfun(@(x) {str2double(x)}, listKOGenes(numIdx));
listKOGenes = [opts.SelectedVariableNames; listKOGenes];
listKOGenes = cellstr(string(listKOGenes));
clear opts

% Retrieve data for measured T-DNAs:
idxKOs = ismember(listKOGenes(:,1), sampleNumbers);
idxKOs(1) = true;
filteredKOGenes = listKOGenes(idxKOs, :);

switch flagFilter
    case 1
        fileName = 'Selected_KO_genes.xlsx';
        % Load genes from model:
        pathModel = fullfile('..', 'IntegrateAraCoreModel', 'OutputFiles', 'OutputModelUniqueBalanced.mat');
        dirModel  = dir(pathModel);
        dirModel = fullfile(dirModel.folder, 'OutputModelUniqueBalanced.mat');
        model = load(dirModel);
        model = model.model;
        modelGenes = model.genes;

        % Now find out how many of the KO genes are included in the model:
        idxKOs = ismember(filteredKOGenes(:,2), modelGenes);
        idxKOs(1) = true;
        selectedKOs = filteredKOGenes(idxKOs, :);
        selectedKOs = [selectedKOs, cell(size(selectedKOs,1),1)];
        selectedKOs{1,5} = 'grRules';
        
        % Find out how many participate as isoenzymes:
        for i = 2:size(selectedKOs,1)
            gene_i = selectedKOs{i,2};
            idxGene_i = contains(model.grRules, gene_i);
            rulesGene_i = unique(model.grRules(idxGene_i));
            selectedKOs{i,5} = rulesGene_i{1};
        end

    case 0
        fileName = 'CompleteSelected_KO_genes.xlsx';
        selectedKOs = filteredKOGenes;
end


if ~exist('flagSave', 'var')
    flagSave = 1;
end

switch flagSave
    case 1
        % Save output variable:
        savePath = fullfile('OutputData', fileName);
        xlswrite(savePath, selectedKOs, 'Selected_KO_genes')
end

end
