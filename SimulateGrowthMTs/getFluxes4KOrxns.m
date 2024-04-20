function fluxes_WtMt = getFluxes4KOrxns(ko_i, fluxesWT, flagSet, flagMethod, rootFolder, indexKO_i)
%##########################################################################
%
% This function uses as input the sampled flux distributions and filters
% the fluxes for the reactions underlying the mutated genes. For the case
% of the ACHR implementation, the fluxes sampled correspond to the closest
% vector that is at steady state, which in turn is closest to a random
% sample from the feasible interval. No further processing steps are
% necessary for the sampled fluxes.
%
% flagSet   - ('wo_relax') without relaxation of lower and upper bounds using lipids sd
%           - ('w_relax') with relaxation of lower and upper bounds using lipids sd
%
% flagMethod - ('distalSampling') Use distal sampling method from Wendering and Nikoloski (https://journals.asm.org/journal/msystems on 12 June 2023 by 95.90.243.53.)
%              ('CHRR') Use Coordinate Hit-and-Run (CHRR) included in the CobraToolbox
%              ('ACHR') Default. Use Artificial Centered Hit-and-Run (ACHR) included in the CobraToolbox  
%
%##########################################################################


switch flagSet
    case {'wo_relax'}
        nameRootFile = '_wo_relax.mat';

    case {'w_relax'}
        nameRootFile = '_w_relax.mat';
end

switch flagMethod
        
    case {'ACHR'}
        % get path for file that contains the samples:
        nameSubfolderModels = 'Models_sampling';
        nameSubfolderSamples = 'Random_samples';
        fileKO_wo = [ko_i, '_modelSampling', nameRootFile];
        dirFileKO_i = dir(fullfile(rootFolder,nameSubfolderModels,fileKO_wo));
        
end


if ~isempty(dirFileKO_i)
    
    switch flagMethod
            
        case {'ACHR'}
            fluxesWT = fluxesWT(:);
            
            % load sampling model for KO_i:
            modelReducedKO_i = load(fullfile(rootFolder, nameSubfolderModels, fileKO_wo));
            modelReducedKO_i = modelReducedKO_i.modelSampling;
            
            % load sampled fluxes for KO_i:
            newFlag = 'w_relax_';
            sampledFluxesKO_i = zeros(numel(modelReducedKO_i.rxns),5000);
            count = 1;
            for j = 1:5
                nameSamplesFile_j = [ko_i, '_ACHR_', newFlag, char(string(j)), '.mat'];
                pathSamples_j = fullfile(rootFolder,nameSubfolderSamples, nameSamplesFile_j);
                samplesMT_j = load(pathSamples_j);
                samplesMT_j = samplesMT_j.sSamples_j;
                sampledFluxesKO_i(:, count:count+(size(samplesMT_j,2)-1)) = samplesMT_j;
                count = count+(size(samplesMT_j,2));
            end
            
            % get the list of reactions catalyzed by the corresponding KO-gene if not provided:
            if nargin < 6
                locusID = char(extractBetween(ko_i, 1, 9));
                indexKO_i = contains(modelReducedKO_i.grRules, locusID);
            end
            
            % get fluxes sampled for selected KO-rxn:
            fluxesKO_i = sampledFluxesKO_i(indexKO_i, :);
            fluxesKO_i = fluxesKO_i(:);
            
    end
    
    % create labels for each flux set:
    labelWT = repmat({'WT'},numel(fluxesWT),1);
    labelKO = repmat({'MT'},numel(fluxesKO_i),1);
    labelsWtMt = cell2table([labelWT;labelKO], "VariableNames",{'Condition'});

    % merge flux samples into one matrix and eliminate rxns not carrying flux:
    fluxes_WtMt = array2table([fluxesWT;fluxesKO_i]);
    fluxes_WtMt.Properties.VariableNames = {'Fluxes'};

    % add the corresponding labels to sets of flux samples:
    fluxes_WtMt = [labelsWtMt,fluxes_WtMt];
    
else
    fluxes_WtMt = [];
end


end

