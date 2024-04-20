function [biomassCoeff, model, LIPIDcoefficients] = getBiomassComponents(accession,day)
%%-------------------------------------------------------------------------
% 'getBiomassComponents' function generates two structures containing the
% stoichiometric coefficients for lipids and other biomass components,
% respectively, using the experimental data provided, for the selected
% accession and day.
% 'DD' refers to the number of days that the plants were subjected to the
% continuous darkness treatment.
%
% USAGE:
%   [biomassCoeff, model, LIPIDcoefficients] = getBiomassComponents(accession,day)
%
% INPUTS:   
%   day       - 0   (Control sample)
%               3   (3DD)
%               6   (6DD)
%               6h  (Col-0, 6h)
%               21h (Col-0, 21h)
%	
%   accession -	1           (when Day = 0)
%               1 up to 282 (when Day = 3)
%               1 up to 284 (when Day = 6)
%               Col0       (when Day = 6h or 21h)
%               All         (when it is necessary to calculate the
%                               stoichiometric coefficients for all the 
%                               accessions in 3DD and 6DD)
%       
% OUTPUTS:
%	model            - COBRA model structure that contains the merged list
%                      of metabolites from the Lipid module and the
%                      ModelTemplate.
%
%   biomassCoeff     - Structure containing the stoichiometric coefficients
%                      for the biomass components: amino acids, soluble 
%                      metabolites, starch, cellulose, nucleotides.
%                      
%   LIPIDcoefficients - Structure containing the stoichiometric
%                       coefficients for the lipid SLIME pseudo-reactions.
%--------------------------------------------------------------------------
%% 1). Get protein-bound amino acids (umol g-1 DW) for samples:
switch day
    case {0, 3}
        dayID = 'D3';
        formatSpec = [repmat('%s', 1, 284), '%[^\n\r]'];

    case 6
        dayID = 'D6';
        formatSpec = [repmat('%s', 1, 286), '%[^\n\r]'];

    case {'6h', '21h'}
        dayID = 'Col0';
        formatSpec = [repmat('%s', 1, 4), '%[^\n\r]'];
end


% Load ultimate amino acids data (umol g-1 DW) for samples:
pathFile = fullfile('SlimeInputData',['UltimateAAs', dayID, '.txt']);
delimiter = {'\t',','};
startRow = 2;
fileID = fopen(pathFile,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);

ultimateAAs = cellstr([dataArray{1:end-1}]);
    
%% 2). Get soluble metabolites (umol g-1 DW) for samples:
pathFile = fullfile('SlimeInputData',['AbsoluteSolMets', dayID, '.txt']); 
delimiter = {'\t',','};
startRow = [2,21];
endRow = [2,30];
fileID = fopen(pathFile,'r');
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end
fclose(fileID);
absoluteSolMets = cellstr([dataArray{1:end-1}]);

% Eliminate empty cells:
emptyAccessions = cellfun(@isempty, absoluteSolMets(2:end,:));
idxEmpty = sum(emptyAccessions) > 1;
absoluteSolMets(:,idxEmpty) = '';

clearvars delimiter startRow endRow formatSpec fileID block dataArrayBlock col dataArray ans;

%% 3). Get nucleotides, starch and cellulose data umol g-1 DW:
pathMets = fullfile('SlimeInputData','Nucleotides_Starch.xlsx');
[~, ~, nucleotidesStarch0_0] = xlsread(pathMets,'Nucleotides','A6:A9');
[~, ~, nucleotidesStarch0_1] = xlsread(pathMets,'Nucleotides','H6:H9');
[~, ~, nucleotidesStarch1_0] = xlsread(pathMets,'Nucleotides','A11:A14');
[~, ~, nucleotidesStarch1_1] = xlsread(pathMets,'Nucleotides','H11:H14');
nucleotides = [nucleotidesStarch0_0,nucleotidesStarch0_1;nucleotidesStarch1_0,nucleotidesStarch1_1;];
nucleotides = string(nucleotides);
nucleotides(ismissing(nucleotides)) = '';
nucleotides = cellstr(nucleotides);

[~, ~, nucleotidesStarch0_0] = xlsread(pathMets,'Starch_Cellulose','A3:F3');
[~, ~, nucleotidesStarch0_1] = xlsread(pathMets,'Starch_Cellulose','A5:F8');
starch_Cellulose = [nucleotidesStarch0_0;nucleotidesStarch0_1];
starch_Cellulose = string(starch_Cellulose);
starch_Cellulose(ismissing(starch_Cellulose)) = '';
starch_Cellulose = cellstr(starch_Cellulose);

clearvars NucleotidesStarch0_0 NucleotidesStarch0_1 NucleotidesStarch1_0 NucleotidesStarch1_1;

%% 4). Get Lipids data:
% Retieve current path and change to directory of 'calculateLipidsCoeff' function:
currentDir = pwd();
pathFun = strrep(currentDir, 'CreateSLIMEmodel', 'ImplementSLIME');
cd(pathFun)

% Obtain name of file that contains lipid coefficients:
switch day
    case 0
        accessionName = num2str(accession);
        dayName = num2str(day);
    case {3, 6}
        accessionName = 'All';
        dayName = num2str(day);
    case {'6h', '21h'}
        accessionName = 'Col0';
        dayName = day;
end

fileName = ['LIPIDcoefficients','_Day',dayName,'_','Accession','-',accessionName,'.mat'];

pathFile = fullfile('..','ImplementSLIMEr','OutputData',fileName);

% Verify if lipid coefficients were already calculated:

if ~exist(pathFile, 'file')
    LIPIDcoefficients = calculateLipidsCoeff(accession,day); 
else
    pathVariable = fullfile('..','ImplementSLIMEr','OutputData',fileName);
    LIPIDcoefficients = load(pathVariable); % Load lipid coefficients into workspace
    LIPIDcoefficients = LIPIDcoefficients.LIPIDcoefficients;
end

cd(currentDir)
clearvars CurrentDir pathFun FindDirFxn DataType fileName pathFile pathVariable

% Retrieve lipids data for selected accession for the data sets 3DD and 6DD:
if isa(accession, 'double') %&& (Day == 3 || Day == 6)
    switch day
        case {3, 6, '6h', '21h'}
            LIPIDcoefficients.Accessions = LIPIDcoefficients.Accessions(accession);
            LIPIDcoefficients.ChainAbundance = LIPIDcoefficients.ChainAbundance(:,accession);
            LIPIDcoefficients.BackboneAbundance = LIPIDcoefficients.BackboneAbundance(:,accession);
            LIPIDcoefficients.TotalLipid_umol = LIPIDcoefficients.TotalLipid_umol(:,accession);
    end
end

%% 5). Eliminate data for accessions that have incomplete data sets:

% First find accessions with complete data sets for lipids and soluble
% metabolites:
[lia,locb] = ismember(LIPIDcoefficients.Accessions,absoluteSolMets(1,:));

LIPIDcoefficients.Accessions        = LIPIDcoefficients.Accessions(lia);
LIPIDcoefficients.ChainAbundance    = LIPIDcoefficients.ChainAbundance(:,lia);
LIPIDcoefficients.BackboneAbundance = LIPIDcoefficients.BackboneAbundance(:,lia);
LIPIDcoefficients.TotalLipid_umol   = LIPIDcoefficients.TotalLipid_umol(lia);

absoluteSolMets = [absoluteSolMets(:,1:2),absoluteSolMets(:,locb((locb > 0)))];

% First find accessions with complete data sets for lipids and amino acids:
[~,locb] = ismember(LIPIDcoefficients.Accessions,ultimateAAs(1,:));
ultimateAAs = [ultimateAAs(:,1:2),ultimateAAs(:,locb(locb > 0))];

%% 6). Load data required to build accession- and condition-specific biomass reactions

% Load list with metabolites required for biomass synthesis:
pathMets = fullfile('SlimeInputData','ListBiomassMets.xlsx');
[~, ~, listBioMets] = xlsread(pathMets,'MetabolitesList','A2:E44');
listBioMets = removeIsNaN(listBioMets);

% Load list with compartments synonyms:
[~, ~, listCompsSyn] = xlsread(pathMets,'CompSynonyms','A2:E4');
listCompsSyn = removeIsNaN(listCompsSyn);

% Load model into workspace:
pathModel = fullfile('SlimeInputData','OutputModelUniqueBalanced.mat');
model = readCbModel(pathModel);

clearvars pathMets idx pathModel FindDir DirPath

% Find compartments abbreviations in model:
comps = unique(cellstr(listBioMets(:,3)));
comps = [comps,cell(size(comps))];

for i = 1:size(comps,1)
    compAbb = comps{i,1};
    compIdx = strcmp(listCompsSyn(:,1),compAbb);
    synonyms = listCompsSyn(compIdx,2:end);
    modelIdx = ismember(model.compNames, synonyms);
    
    if sum(modelIdx) == 1
        GetAbb = model.comps{modelIdx};
        comps{i,2} = ['[',GetAbb,']'];
    end   
end

% Find metabolites IDs used in model:
listBioMets = cellstr([listBioMets,strings(size(listBioMets,1),1)]);

for i = 1:size(listBioMets,1)
    compAbb = listBioMets{i,3}; % Get compartment ID
    compIdx = strcmp(comps(:,1),compAbb);
    getCompAbb = comps{compIdx,2};
    getKEGGID = listBioMets{i,5}; % Get KEGG ID
    getMetName = listBioMets{i,1}; % Get Metabolite name
       
    if isempty(getKEGGID)
        modelMetID = contains(model.mets,getCompAbb).*contains(model.mets,getMetName)...
            .*~contains(model.mets,{'chain','backbone'})== 1;
    else
        modelMetID = strcmp(model.metKEGGID, getKEGGID).*contains(model.mets,getCompAbb) == 1;
        if sum(modelMetID) > 1
            modelMetID = strcmp(model.metKEGGID, getKEGGID).*contains(model.mets,getCompAbb)...
                .*contains(model.mets,getMetName) == 1;
        end
    end
    
    listBioMets{i,6} = model.mets(modelMetID);  
end

clearvars i idx pathModel FindDir DirPath CompAbb CompIdx Synonyms...
    ModelIdx GetAbb GetCompAbb ModelMetID GetKEGGID GetMetName

%% 7). Build matrix with stoichiometric coefficients of biomass reactions:

% Create array structure and add coefficients data:
biomassCoeff.Mets = cellstr(string(listBioMets(:,6))); % Create metabolite names list
biomassCoeff.S = zeros(size(biomassCoeff.Mets,1),size(LIPIDcoefficients.Accessions,2));

if day == 0
    starchIdx = 2;
elseif day == 3
    starchIdx = 3;
elseif day == 6
    starchIdx = 4;
elseif strcmp(day,'6h')
    starchIdx = 5;
elseif strcmp(day,'21h')
    starchIdx = 6;
end
        
for i = 1:size(listBioMets,1)
    metName = listBioMets{i,2};
    metCategory = listBioMets{i,4};
    
    for j = 1:size(LIPIDcoefficients.Accessions,2)
        nameAccession = LIPIDcoefficients.Accessions{1,j};
        
        switch metCategory
            case {'AA'}
                switch day
                    case 0
                        findAcession = 2;
                    case {3, 6, '6h', '21h'}
                        findAcession = strcmp(ultimateAAs(1,:), nameAccession);
                end
                aaIdx = strcmp(ultimateAAs(:,1),metName); % Find accesion data
                abundance = str2double(string(ultimateAAs(aaIdx,findAcession)));
                biomassCoeff.S(i,j) = abundance;
                
            case {'Nucleotides'}
                nucIdx = strcmp(nucleotides(:,1),metName); % Find accesion data
                abundance = str2double(string(nucleotides(nucIdx,2)));
                biomassCoeff.S(i,j) = abundance;
                
            case {'SolMets'}
                switch day
                    case 0
                        findAcession = 2;
                    case {3, 6, '6h', '21h'}
                        findAcession = strcmp(absoluteSolMets(1,:), nameAccession);
                end
                solMetIdx = strcmp(absoluteSolMets(:,1),metName); % Find accesion data
                abundance = str2double(string(absoluteSolMets(solMetIdx,findAcession)));
                biomassCoeff.S(i,j) = abundance;
                
            case {'Cel-Starch-Others'}
                othersIdx = strcmp(starch_Cellulose(:,1),metName); % Find accesion data
                abundance = str2double(string(starch_Cellulose(othersIdx,starchIdx)));
                biomassCoeff.S(i,j) = abundance;
                
            case {'Lipid'}
                biomassCoeff.S(i,j) = LIPIDcoefficients.TotalLipid_umol(1,j);
        end
    end
end

biomassCoeff.S(:,all(~any(biomassCoeff.S),1)) = []; % Remove columns with zero entries


clearvars Abundance i j IdxEmpty MetID MetIdx NucIdx OthersIdx...
    SolMetsIdx

end
