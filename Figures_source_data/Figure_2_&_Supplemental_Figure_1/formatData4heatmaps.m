function [prunedFC_mpi,prunedFC_lusk,prunedScaffold_mpi,prunedScaffold_lusk] = formatData4heatmaps
%##########################################################################
%
% This script is used to process the lipid dataset to create a heatmap
% where lipid species of particular interest are highlighted: (i) species
% that show changes acccording to literature descriptions, (ii) species
% showing significant changes in relative abundances).
%
%##########################################################################


% ###
% ##### 1). Process list of mutants dataset from MPI:
% ###
%pathDatasetConcordant = 'T_test_lipidome_concordantMTs_16-01-2024_edited.xlsx';
%pathDatasetOutlier = 'T_test_lipidome_outlierMTs_16-01-2024_edited.xlsx';

pathDatasetConcordant = 'T_test_lipidome_concordantMTs.xlsx';
pathDatasetOutlier = 'T_test_lipidome_outliersMTs.xlsx';
pathDatasetUnknown = 'T_test_lipidome_unknownMTs.xlsx';


lipidIDs_mpi = readcell('List_locus_and_lipids.xlsx', 'Sheet', 'lipidIDs_mpi');
mpiConcordant = {'AT1G08510','AT1G74960','AT2G32260','AT2G43710','AT3G18000','AT5G46290'};
mpiOutlier = {'AT1G73600','AT1G74320','AT1G77590','AT3G02610','AT3G07690','AT3G15850','AT4G00550','AT4G22340','AT4G25970','AT4G27030','AT4G33030','AT5G15530'};

colNames = [repmat({'concordant'}, 1, numel(mpiConcordant)),repmat({'outlier'}, 1, numel(mpiOutlier))];
locus_mpi = [mpiConcordant, mpiOutlier];

scaffoldMPI = zeros(numel(lipidIDs_mpi),numel(colNames));
scaffoldMPI = createSignificanceMatrix(mpiConcordant,mpiOutlier,scaffoldMPI, pathDatasetConcordant, pathDatasetOutlier,lipidIDs_mpi, []);

lipidIDs_mpi = [{'mutant class';'species'};lipidIDs_mpi];
scaffoldMPI = [locus_mpi;colNames;cellstr(string(scaffoldMPI))];
scaffoldMPI = [lipidIDs_mpi,scaffoldMPI];




% ###
% ##### 2). Process list of MTs from Lusk:
% ###
lipidIDs_lusk = readcell('List_locus_and_lipids.xlsx', 'Sheet', 'lipidIDsLusk');
luskConcordant = {'AT1G15080','AT1G32200','AT1G80950','AT3G14360','AT3G15730','AT3G26840','AT3G57140','AT5G04040','AT5G57190'};
luskOutlier = {'AT1G02390','AT1G51260','AT1G75020','AT2G19450','AT3G18850','AT5G42870','AT5G66450'};

colNames = [repmat({'concordant'}, 1, numel(luskConcordant)),repmat({'outlier'}, 1, numel(luskOutlier))];
locus_lusk = [luskConcordant, luskOutlier];

scaffoldLusk = zeros(numel(lipidIDs_lusk),numel(colNames));
scaffoldLusk = createSignificanceMatrix(luskConcordant,luskOutlier,scaffoldLusk, pathDatasetConcordant, pathDatasetOutlier,lipidIDs_lusk, []);

lipidIDs_lusk = [{'mutant class';'species'};lipidIDs_lusk];
scaffoldLusk = [locus_lusk;colNames;cellstr(string(scaffoldLusk))];
scaffoldLusk = [lipidIDs_lusk,scaffoldLusk];



% ###
% ##### 3). Create list of locus and lipids to include in heatmaps:
% ###
% load fold-change computed for lipids measured for MTs:
currentFolder = pwd;
dirFx = strrep(currentFolder,'Figures_source_data\Figure_2_&_Supplemental_Figure_1', 'SimulateGrowthMTs');
oldDir = cd(dirFx);
[~, ~, ~, ~, processedData_mpi] = constrainUnmeasuredSpecies('mpi', 0, 0);
[~, ~, ~, ~, processedData_lusk] = constrainUnmeasuredSpecies('lusk', 0, 0);
cd(oldDir)
fc_mpi = processedData_mpi.processedProfiles.FCaverages;
fc_mpi(2:end,1) = extractAfter(fc_mpi(2:end,1), '] ');
fc_mpi(:,1) = strrep(fc_mpi(:,1), 'PG 34:4 ', 'PG 34:4');%adjust this one manually
fc_lusk = processedData_lusk.processedProfiles.FCaverages;
fc_lusk(2:end,1) = extractAfter(fc_lusk(2:end,1), '] ');


% create list of locus to keep:
listSelectedLocus_mpi = [mpiConcordant,mpiOutlier];

listSelectedLocus_lusk = [luskConcordant,luskOutlier];


% eliminate locus that won't be considered in the heatmaps:
prunedFC_mpi = pruneDataset(fc_mpi, listSelectedLocus_mpi, 'locus', 0);
prunedScaffold_mpi = pruneDataset(scaffoldMPI, listSelectedLocus_mpi, 'locus', 1);

prunedFC_lusk = pruneDataset(fc_lusk, listSelectedLocus_lusk, 'locus', 0);
prunedScaffold_lusk = pruneDataset(scaffoldLusk, listSelectedLocus_lusk, 'locus', 1);


% load list of lipids to keep:
pathLocus = 'List_locus_and_lipids.xlsx';
lipids2keep_mpi = cellstr(string(readcell(pathLocus, 'Sheet', 'list_mpi')));
lipids2keep_lusk = cellstr(string(readcell(pathLocus, 'Sheet', 'list_lusk')));

% eliminate lipid species that won't be considered in the heatmaps:
prunedFC_mpi = pruneDataset(prunedFC_mpi, lipids2keep_mpi, 'species', 0);
prunedScaffold_mpi = pruneDataset(prunedScaffold_mpi, lipids2keep_mpi, 'species', 1);

prunedFC_lusk = pruneDataset(prunedFC_lusk, lipids2keep_lusk, 'species', 0);
prunedScaffold_lusk = pruneDataset(prunedScaffold_lusk, lipids2keep_lusk, 'species', 1);

% carry out manual adjustments to Lusk dataset:
idxPE = strcmp(prunedFC_lusk(:,1), 'PE 42:3');
avPE = mean(str2double(string(prunedFC_lusk(idxPE,2:end))));
idxPE = find(idxPE);
prunedFC_lusk(idxPE(1),2:end) = cellstr(string(avPE));
prunedFC_lusk(idxPE(2),:) = '';

idxPE = find(strcmp(prunedScaffold_lusk(:,1), 'PE 42:3'));
prunedScaffold_lusk(idxPE(2),:) = '';
idxPI = find(strcmp(prunedScaffold_lusk(:,1), 'PI 34:3'));
prunedScaffold_lusk(idxPI(2),:) = '';


% create tables and save in .xlsx file:
path2save = 'Lipid_profiles4heatmaps.xlsx';
colsMPI = prunedFC_mpi(1,:);
rowsMPI = cell2table(prunedFC_mpi(2:end,1));
rowsMPI.Properties.VariableNames = "locus";
fcsMPI = array2table(str2double(string(prunedFC_mpi(2:end,2:end))));
fcsMPI = [rowsMPI,fcsMPI];
fcsMPI.Properties.VariableNames = colsMPI;
writetable(fcsMPI,path2save,'Sheet', 'mpi_dataset')

ps_mpi = array2table(str2double(string(prunedScaffold_mpi(3:end,2:end))));
ps_mpi = [rowsMPI,ps_mpi];
ps_mpi.Properties.VariableNames = colsMPI;
writetable(ps_mpi,path2save,'Sheet', 'scaffold_asterisk_mpi')


colsLusk = prunedFC_lusk(1,:);
rowsLusk = cell2table(prunedFC_lusk(2:end,1));
rowsLusk.Properties.VariableNames = "locus";
fcsLusk = array2table(str2double(string(prunedFC_lusk(2:end,2:end))));
fcsLusk = [rowsLusk,fcsLusk];
fcsLusk.Properties.VariableNames = colsLusk;
writetable(fcsLusk,path2save,'Sheet', 'lusk_dataset')

ps_lusk = array2table(str2double(string(prunedScaffold_lusk(3:end,2:end))));
ps_lusk = [rowsLusk,ps_lusk];
ps_lusk.Properties.VariableNames = colsLusk;
writetable(ps_lusk,path2save,'Sheet', 'scaffold_asterisk_lusk')




% ###
% ##### 4). Create list of locus and lipids that were not included in the main heatmap:
% ###
% This below set will be used to create a supplementary figure!

%            'locisID'   'data set'  'details'
listLocus = {'AT3G48780'    'mpi'    'not meassured';
            'AT5G23670'    'mpi'    'not meassured';
            'AT1G06520'    'mpi'    'not obvious';
            'AT5G19200'    'mpi'    'not meassured';
            'AT4G23850'    'mpi'    'PGPS1';
            'AT1G64670'    'lusk'    'not meassured';
            'AT2G37940'    'lusk'    'not meassured';
            'AT3G05630'    'lusk'    'mutant is not characterized';
            'AT3G56940'    'lusk'    'not meassured';
            'AT4G26910'    'lusk'    'not meassured';
            'AT5G25370'    'lusk'    'differences among lipidome measured and phenotype'; % flux of this MT agreed with pheno, but changes in lipidome don´t agree with phenotipe reported
            'AT5G37300'    'lusk'    'not meassured';
            'AT1G01610'    'lusk'    'not meassured';
            'AT2G39290'    'lusk'    'lipidome agree / artifact sampling';
            'AT2G42010'    'lusk'    'scarce phenotypic data available';
            'AT4G00400'    'lusk'    'not meassured';
            'AT4G18550'    'lusk'    'not obvious';
            'AT5G03080'    'lusk'    'not obvious';
            'AT1G72520'    'lusk'    'LOX4';
            'AT2G38110'    'lusk'    'GPAT6';
            'AT3G03540'    'lusk'    'unknown';
            'AT3G57650'    'lusk'    'LPAT2'};

        
% get lipid profiles for selected sets of mutants:
heatmap2mpi = cell(size(fc_mpi,1),5);
heatmap2mpi(:,1) = fc_mpi(:,1);
heatmap2lusk = cell(size(fc_lusk,1),14);
heatmap2lusk(:,1) = fc_lusk(:,1);
count_mpi = 2;
count_lusk = 2;

for i = 1:size(listLocus)
    locus_i = listLocus{i,1};
    dataSet = listLocus{i,2};

    switch dataSet
        case {'mpi'}
            idxLocus = strcmp(fc_mpi(1,:), locus_i);
            heatmap2mpi(:,count_mpi) = fc_mpi(:,idxLocus);
            count_mpi = count_mpi+1;
            
        case {'lusk'}
            idxLocus = strcmp(fc_lusk(1,:), locus_i);
            heatmap2lusk(:,count_lusk) = fc_lusk(:,idxLocus);
            count_lusk = count_lusk+1;
    end
    
end

% perform manual edits to make easier constructing heatmaps in R:
% a). mpi data set: delete species with low confidence or uninteresting species.
species2delete = {'DGDG 32:', 'MGDG 32:', 'MGMG '};

for i = 1:numel(species2delete)
    species_i = species2delete{i};
    idxSpecies_i = contains(heatmap2mpi(:,1), species_i);
    heatmap2mpi(idxSpecies_i,:) = '';
end

% create scaffold to assign significance asteriks for mpi data set:
lipidIDs_mpi2 = readcell('List_locus_and_lipids.xlsx', 'Sheet', 'lipidIDs_mpi_HM2');
mpiConcordant2 = listLocus(1:2,1);
mpiOutlier2 = listLocus(3:5,1);
colNames = [repmat({'concordant'}, 1, numel(mpiConcordant2)),repmat({'outlier'}, 1, numel(mpiOutlier2))];
locus_mpi2 = [mpiConcordant2', mpiOutlier2'];

scaffoldMPI2 = zeros(numel(lipidIDs_mpi2),numel(colNames));
scaffoldMPI2 = createSignificanceMatrix(mpiConcordant2,mpiOutlier2,scaffoldMPI2,pathDatasetConcordant,pathDatasetOutlier,lipidIDs_mpi2,pathDatasetUnknown);
lipidIDs_mpi2 = [{'mutant class';'species'};lipidIDs_mpi2];
scaffoldMPI2 = [locus_mpi2;colNames;cellstr(string(scaffoldMPI2))];
scaffoldMPI2 = [lipidIDs_mpi2,scaffoldMPI2];

% create tables and save results for building supplementary mpi heatmap:
colsMPI2 = heatmap2mpi(1,:);
rowsMPI2 = cell2table(heatmap2mpi(2:end,1));
rowsMPI2.Properties.VariableNames = "locus";
fcsMPI2 = array2table(str2double(string(heatmap2mpi(2:end,2:end))));
fcsMPI2 = [rowsMPI2,fcsMPI2];
fcsMPI2.Properties.VariableNames = colsMPI2;
writetable(fcsMPI2,path2save,'Sheet', 'second_heatmap_mpi')

ps2_mpi = array2table(str2double(string(scaffoldMPI2(3:end,2:end))));
ps2_mpi = [rowsMPI2,ps2_mpi];
ps2_mpi.Properties.VariableNames = colsMPI2;
writetable(ps2_mpi,path2save,'Sheet', 'second_scaffold_mpi')


% b). lusk data set: compute average for duplicated species.
idxPE = strcmp(heatmap2lusk(:,1), 'PE 42:3');
avPE = mean(str2double(string(heatmap2lusk(idxPE,2:end))));
idxPE = find(idxPE);
heatmap2lusk(idxPE(1),2:end) = cellstr(string(avPE));
heatmap2lusk(idxPE(2),:) = '';

% create scaffold to assign significance asteriks for lusk data set:
lipidIDs_lusk2 = readcell('List_locus_and_lipids.xlsx', 'Sheet', 'lipidIDs_lusk_HM2');
luskConcordant2 = listLocus(6:12,1);
luskOutlier2 = listLocus(13:end,1);
colNames = [repmat({'concordant'}, 1, numel(luskConcordant2)),repmat({'outlier'}, 1, numel(luskOutlier2))];
locus_lusk2 = [luskConcordant2', luskOutlier2'];

scaffoldLusk2 = zeros(numel(lipidIDs_lusk2),numel(colNames));
scaffoldLusk2 = createSignificanceMatrix(luskConcordant2,luskOutlier2,scaffoldLusk2,pathDatasetConcordant,pathDatasetOutlier,lipidIDs_lusk2,pathDatasetUnknown);
lipidIDs_lusk2 = [{'mutant class';'species'};lipidIDs_lusk2];
scaffoldLusk2 = [locus_lusk2;colNames;cellstr(string(scaffoldLusk2))];
scaffoldLusk2 = [lipidIDs_lusk2,scaffoldLusk2];


% create tables and save results for building supplementary mpi heatmap:
colsLusk2 = heatmap2lusk(1,:);
rowsLusk2 = cell2table(heatmap2lusk(2:end,1));
rowsLusk2.Properties.VariableNames = "locus";
fcsLusk2 = array2table(str2double(string(heatmap2lusk(2:end,2:end))));
fcsLusk2 = [rowsLusk2,fcsLusk2];
fcsLusk2.Properties.VariableNames = colsLusk2;
writetable(fcsLusk2,path2save,'Sheet', 'second_heatmap_lusk')

ps2_lusk = array2table(str2double(string(scaffoldLusk2(3:end,2:end))));
ps2_lusk = [rowsLusk2,ps2_lusk];
ps2_lusk.Properties.VariableNames = colsLusk2;
writetable(ps2_lusk,path2save,'Sheet', 'second_scaffold_lusk')


end

%%
function prunedData = pruneDataset(fc, list2filter, flagCase, flagScaffold)

switch flagCase
    case {'locus'}
        prunedData = cell(size(fc,1), numel(list2filter));
        for i = 1:numel(list2filter)
            locus_i = list2filter{i};
            idxLocus_i = strcmp(fc(1,:), locus_i);
            prunedData(:,i) = fc(:,idxLocus_i);
        end
        prunedData = [fc(:,1),prunedData];
       

    case {'species'}
        idxLipids = ismember(fc(:,1), list2filter);
        if flagScaffold == 1
            idxLipids(1:2) = true;
        else
            idxLipids(1) = true;
        end
        prunedData = fc(idxLipids,:);
        
end

end



%%
function scaffoldMPI = createSignificanceMatrix(mpiConcordant,mpiOutlier,scaffoldMPI, pathDatasetConcordant, pathDatasetOutlier,lipidIDs_mpi, pathDatasetUnknown)

%luskUnknown = {'AT3G03540','AT1G72520','AT2G38110','AT3G57650'};
%mpiUnknown = {'AT4G23850'};  

count = 1;
% check profiles for concordant MTs:
for i = 1:numel(mpiConcordant)
    locus_i = mpiConcordant{i};
    
    switch locus_i
        case {'AT3G03540','AT1G72520','AT2G38110','AT3G57650','AT4G23850'}
            profileLocus_i = cellstr(string(readcell(pathDatasetUnknown, 'Sheet', locus_i)));
        otherwise
            profileLocus_i = cellstr(string(readcell(pathDatasetConcordant, 'Sheet', locus_i)));
    end
    
    profileLocus_i(2:end,1) = extractAfter(profileLocus_i(2:end,1), '] ');
    idxHits = strcmp(profileLocus_i(:,12), 'Y');
    listHits = (profileLocus_i(idxHits, 1:2));
    
    % match to scaffold list lipids that contribute to phenotype:
    idxMatchedPheno = ismember(listHits(:,1),lipidIDs_mpi).*strcmp(listHits(:,2),'0');
    idxHitsScaffold = ismember(lipidIDs_mpi, listHits(idxMatchedPheno==1,1));
    scaffoldMPI(idxHitsScaffold,count) = 1;
    
    % match to scaffold lipids that contribute to pheno and changed significantly:
    idxMatchedPheno = ismember(listHits(:,1),lipidIDs_mpi).*strcmp(listHits(:,2),'1');
    idxHitsScaffold = ismember(lipidIDs_mpi, listHits(idxMatchedPheno==1,1));
    scaffoldMPI(idxHitsScaffold,count) = 2;
    count = count+1;
end

% check profiles for outliers MTs:
for i = 1:numel(mpiOutlier)
    locus_i = mpiOutlier{i};
    
    switch locus_i
        case {'AT3G03540','AT1G72520','AT2G38110','AT3G57650','AT4G23850'}
            profileLocus_i = cellstr(string(readcell(pathDatasetUnknown, 'Sheet', locus_i)));
        otherwise
            profileLocus_i = cellstr(string(readcell(pathDatasetOutlier, 'Sheet', locus_i)));
    end
    
    profileLocus_i(2:end,1) = extractAfter(profileLocus_i(2:end,1), '] ');
    idxHits = strcmp(profileLocus_i(:,12), 'Y');
    listHits = (profileLocus_i(idxHits, 1:2));
    
    % match to scaffold list lipids that contribute to phenotype:
    idxMatchedPheno = ismember(listHits(:,1),lipidIDs_mpi).*strcmp(listHits(:,2),'0');
    idxHitsScaffold = ismember(lipidIDs_mpi, listHits(idxMatchedPheno==1,1));
    scaffoldMPI(idxHitsScaffold,count) = 1;
    
    % match to scaffold lipids that contribute to pheno and changed significantly:
    idxMatchedPheno = ismember(listHits(:,1),lipidIDs_mpi).*strcmp(listHits(:,2),'1');
    idxHitsScaffold = ismember(lipidIDs_mpi, listHits(idxMatchedPheno==1,1));
    scaffoldMPI(idxHitsScaffold,count) = 2;
    count = count+1;
end


end