function [feasibleModels,unfeasibleModels] = sampleDist4mutants(growthMTs, modelsMT, flagMethod, flagModelRed)
%sampleDist4mutants(growthMTs, modelsMT, flagModels, flagMethod) % oldFx
%##########################################################################
%
% USAGE: 
%   [fva_MTs] = sampleDist4mutants(growthMTs, modelsMT, flagModels, flagMethod)
%
% INPUTS:
%   flagMethod  -  ('ACHR') Default
%
%   flagModelRed - (false) Default. Sampling is carried out without reducing the model
%                  (true) Reduce model by eliminating reactions not carrying flux
%
%##########################################################################


% identify mutants for which the constrained model was unfeasible and
% split accordingly:
idxFeasible = ~cellfun(@isempty, growthMTs(:,2));
feasibleModels = growthMTs(idxFeasible, :);
unfeasibleModels = growthMTs(idxFeasible==0, :);
dirFVAresults = dir(fullfile('..', 'SamplingResults', 'FVA_results', 'tempVariable.mat'));
pathFVAresults = dirFVAresults.folder;
flagRelaxation = '_w_relax.mat';
dirSamples = strrep(pathFVAresults, 'FVA_results', 'Random_samples');
dirSamplingFiles = fullfile(strrep(pathFVAresults, 'FVA_results', 'FilesSamplingProcedure'), 'W_relaxation');
dirSamplingModels = strrep(pathFVAresults, 'FVA_results', 'Models_sampling');

if nargin < 4
    flagModelRed = false;
end

% get feasible models, perform FVA if requested and then do the sampling:
for i = 2:size(feasibleModels,1)
    % get model built for mutant:
    mutant_i = feasibleModels{i,1};
    modelMT_i = modelsMT.(mutant_i);
    
    switch flagModelRed
        case 0 % perform FVA
            % check if fva analysis was done:
            nameFile = ['fva_', mutant_i, flagRelaxation];
            bioOptMT_i = str2double(string(feasibleModels(i,2)));
            if isempty(dir(fullfile(pathFVAresults,nameFile)))
                % constrain biomass reaction with corresponding optimum biomass:
                modelMT_i.ub(modelMT_i.c~=0) = bioOptMT_i;
                [minFluxMT_i, maxFluxMT_i] = fluxVariability(modelMT_i,100,'max');
                fva_MT.minFlux = minFluxMT_i;
                fva_MT.maxFlux = maxFluxMT_i;
                save(fullfile(pathFVAresults,nameFile), 'fva_MT')

            else
                fva_MT = load(fullfile(pathFVAresults,nameFile));
                fva_MT = fva_MT.fva_MT;
                minFluxMT_i = fva_MT.minFlux;
                maxFluxMT_i = fva_MT.maxFlux;
            end
    
        case 1 % do not perform FVA since this will be done during model reducction
            minFluxMT_i = '';
            maxFluxMT_i = '';
    end
    
    % select method for sampling the solution space:
    modelMT_i.ub(modelMT_i.c~=0) = 1000;
    bioOptMT_i = str2double(string(feasibleModels(i,2)));

    switch flagMethod

        case {'ACHR'}
            % create name for file to save sampling results:
            newFlag = strrep(flagRelaxation, '.mat', '_1.mat');
            nameExistingFile = [mutant_i, '_ACHR', newFlag];
            pathFileSamples = fullfile(dirSamples, nameExistingFile);
            fileName = extractBefore(nameExistingFile, '_1.mat');

            % run 'sampleCbModel' if the analysis was not carried out:
            if isempty(dir(pathFileSamples))
                options.nStepsPerPoint = 200;
                options.nPointsReturned = 5000;
                options.nWarmupPoints = 15000; %previous 15000
                options.nFiles = 10;
                options.nPointsPerFile = 1000;
                options.nFilesSkipped = 2;
                options.minFlux = minFluxMT_i;
                options.maxFlux = maxFluxMT_i;
                options.bioOpt = bioOptMT_i;
                options.mutant_i = mutant_i;
                options.flagRelaxation = flagRelaxation;
                options.flagModelReduction = flagModelRed;
                options.routeFolder = dirSamplingFiles;
                [modelSampling, samplesMT_i, ~] = sampleCbModel_modified(modelMT_i, fileName, 'ACHR', options);
                %[modelSampling,samples,volume] = sampleCbModel_modified(model, sampleFile, samplerName, options, modelSampling)
                
                %split samples to reduce file size and save samples:
                count = 1;
                for j = 1:5
                    sSamples_j = samplesMT_i(:,count:(1000*j));
                    nameSamples_j = [fileName, '_', char(string(j)), '.mat'];
                    pathSamples_j = fullfile(dirSamples, nameSamples_j);
                    save(pathSamples_j, 'sSamples_j')
                    count = count+(1000);
                end

                % save the model sampling for later use!
                nameModelSampling = [mutant_i, '_modelSampling', flagRelaxation];
                pathSamplingModels = fullfile(dirSamplingModels,nameModelSampling);
                save(pathSamplingModels, 'modelSampling')
                
            end
    end
end

end