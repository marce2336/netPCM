function [statsLipidsOutliers,statsLipidsConcordants,statsLipidsUnknown] = analyzeLipidomes
%##########################################################################
%
% Analyze the lipid profiles measured for each of the mutants (MTs) to
% identify case-by-case the change in the direction of the abundances for
% the lipid species measured compared to the WT-plants
%
%##########################################################################

% ###
% ##### 1). Load lipid profiles for outliers identified for each dataset: 
% ###
dirFx = pwd;
dirCusFx = strrep(dirFx, 'Figures_source_data\Figure_2_&_Supplemental_Figure_1', 'SimulateGrowthMTs');
oldDir = cd(dirCusFx);
[~, ~, ~, ~, processedData_mpi] = constrainUnmeasuredSpecies('mpi', 0, 0);
[~, ~, ~, ~, processedData_lusk] = constrainUnmeasuredSpecies('lusk', 0, 0);
cd(oldDir)




% ###
% ##### 2). Analyze lipid profiles of outliers with published data: 
% ###

% list of outliers with published lipid profiles:
outliers_w_lipids = {'AT1G06520' 'mpi' 'all';
    'AT1G73600' 'mpi' 'all';
    'AT1G74320' 'mpi' 'all';
    'AT1G77590' 'mpi' 'all';
    'AT3G02610' 'mpi' 'all';
    'AT3G07690' 'mpi' 'all';
    'AT3G15850' 'mpi' 'MGDG,DGDG';
    'AT4G00550' 'mpi' 'MGDG,DGDG';
    'AT4G22340' 'mpi' 'all';
    'AT4G25970' 'mpi' 'PE';
    'AT4G27030' 'mpi' 'all';
    'AT4G33030' 'mpi' 'SQDG';
    'AT5G15530' 'mpi' 'all';
    'AT5G19200' 'mpi' 'all';
    'AT1G01610' 'lusk' 'all';
    'AT1G02390' 'lusk' 'all';
    'AT1G51260' 'lusk' 'all';
    'AT1G75020' 'lusk' 'all';
    'AT2G19450' 'lusk' 'all';
    'AT2G39290' 'lusk' 'all';
    'AT2G42010' 'lusk' 'all';
    'AT3G18850' 'lusk' 'all';
    'AT4G00400' 'lusk' 'all';
    'AT4G18550' 'lusk' 'all';
    'AT5G03080' 'lusk' 'all';
    'AT5G42870' 'lusk' 'all';
    'AT5G66450' 'lusk' 'all'};

pathSave = fullfile(dirFx,'T_test_lipidome_outliersMTs.xlsx');
statsLipidsOutliers = compareLipidProfiles(outliers_w_lipids,processedData_mpi,processedData_lusk,pathSave);





% ###
% ##### 3). Analyze lipid profiles of concordant mutants: 
% ###

concordant_w_lipids = {'AT1G08510'	'mpi'	'DG,DGDG,MGDG,SQDG,PI,PG,PE,PC,TG';
    'AT1G74960'	'mpi'	'DG,DGDG,MGDG,SQDG,PI,PG,PE,PC';
    'AT2G32260'	'mpi'	'PC,PS,MGDG';
    'AT2G43710'	'mpi'	'DG,DGDG,MGDG,SQDG,PI,PG,PE,PC,TG';
    'AT3G18000'	'mpi'	'PC,LPC';
    'AT5G23670'	'mpi'	'all';
    'AT3G48780'	'mpi'	'all';
    'AT5G23670'	'mpi'	'all';
    'AT5G46290'	'mpi'	'all';
    'AT1G15080' 'lusk' 'PA,DG'
    'AT1G32200'	'lusk'	'LPA,PA,DG,PG,MGDG,DGDG';
    'AT1G64670'	'lusk'	'all';
    'AT1G80950'	'lusk'	'all';
    'AT2G37940'	'lusk'	'all';
    'AT3G05630'	'lusk'	'all';
    'AT3G14360'	'lusk'	'all';
    'AT3G15730'	'lusk'	'PA,PG';
    'AT3G26840'	'lusk'	'all';
    'AT3G56940'	'lusk'	'all';
    'AT3G57140'	'lusk'	'all';
    'AT4G26910'	'lusk'	'all';
    'AT5G04040'	'lusk'	'all';
    'AT5G25370'	'lusk'	'all';
    'AT5G37300'	'lusk'	'all';
    'AT5G57190'	'lusk'	'all';
    };

pathSave2 = fullfile(dirFx,'T_test_lipidome_concordantMTs.xlsx');
statsLipidsConcordants = compareLipidProfiles(concordant_w_lipids,processedData_mpi,processedData_lusk,pathSave2);




% ###
% ##### 4). Analyze lipid profiles of mutants lacking phenotypic data: 
% ###

unknown_w_lipids = {'AT3G03540' 'lusk' 'all';
    'AT4G23850' 'mpi' 'all';
    'AT1G72520' 'lusk' 'all';
    'AT2G38110' 'lusk' 'all';
    'AT3G57650' 'lusk' 'all'};
 
pathSave3 = fullfile(dirFx,'T_test_lipidome_unknownMTs.xlsx');
statsLipidsUnknown = compareLipidProfiles(unknown_w_lipids,processedData_mpi,processedData_lusk,pathSave3);



% ###
% ##### 5). Identify active reactions: 
% ###
oldDir = cd(dirCusFx);
referenceFluxes_mpi = processedData_mpi.referenceFluxes;
referenceFluxes_mpi = sumIncomingFluxes(referenceFluxes_mpi, 0);
lipidIDs_mpi = processedData_mpi.processedProfiles.bioRepsMTs.lipidIDs(2:end,1);
lipidIDs_mpi = extractAfter(lipidIDs_mpi, '] ');
lipidIDs_mpi = strrep(lipidIDs_mpi, ' ', '_');
lipidIDs_mpi = strrep(lipidIDs_mpi, ':', '_');
lipidIDs_mpi = [lipidIDs_mpi, repmat({'y'}, numel(lipidIDs_mpi),1)];
idxModelIDs = cellfun(@isempty, processedData_mpi.processedProfiles.bioRepsMTs.lipidIDs(2:end,2));
lipidIDs_mpi(idxModelIDs,2) = {'n'};
idxHits = find(idxModelIDs==0);
lipidIDs_mpi = [lipidIDs_mpi, repmat({'n.a.'}, size(lipidIDs_mpi,1),1)];
for i = 1:numel(idxHits)
    species_i = lipidIDs_mpi{idxHits(i)};
    fsr_i = referenceFluxes_mpi.(species_i).fsr;
    lipidIDs_mpi(idxHits(i),3) = cellstr(string(fsr_i));
end


referenceFluxes_lusk = processedData_lusk.referenceFluxes;
referenceFluxes_lusk = sumIncomingFluxes(referenceFluxes_lusk, 0);
lipidIDs_lusk = processedData_lusk.processedProfiles.bioRepsMTs.lipidIDs(2:end,2);
lipidIDs_lusk = extractAfter(lipidIDs_lusk, '] ');
lipidIDs_lusk = strrep(lipidIDs_lusk, ' ', '_');
lipidIDs_lusk = strrep(lipidIDs_lusk, ':', '_');
lipidIDs_lusk = [lipidIDs_lusk, repmat({'y'}, numel(lipidIDs_lusk),1)];
idxModelIDs = cellfun(@isempty, processedData_lusk.processedProfiles.bioRepsMTs.lipidIDs(2:end,3));
lipidIDs_lusk(idxModelIDs,2) = {'n'};
idxHits = find(idxModelIDs==0);
lipidIDs_lusk = [lipidIDs_lusk, repmat({'n.a.'}, size(lipidIDs_lusk,1),1)];
for i = 1:numel(idxHits)
    species_i = lipidIDs_lusk{idxHits(i)};
    fsr_i = referenceFluxes_lusk.(species_i).fsr;
    lipidIDs_lusk(idxHits(i),3) = cellstr(string(fsr_i));
end

cd(oldDir)


end