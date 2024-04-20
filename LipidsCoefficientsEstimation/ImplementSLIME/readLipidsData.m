function data = readLipidsData(Accession,Day)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% modified from function: data = readEjsingData(i), created by Benjamín J.
% Sánchez (Last update: 2018-03-20).
%
% Usage
%
%       data = readLipidsData(Accession,Day)
%
% Inputs
%
%       Accession - 1           (when Day = 0)
%                   1 up to 282 (when Day = 3)
%                   1 up to 284 (when Day = 6)
%                   Col0       (when Day = 6h or 21h)
%                   All         (when it is necessary to calculate the
%                               stoichiometric coefficients for all the 
%                               accessions in 3DD and 6DD)
%       
%       Day -   0   (Control sample)
%               3   (3DD)
%               6   (6DD)
%               6h  (Col-0, 6h)
%               21h (Col-0, 21h)
%               Herein 'DD' refers to the number of days that the plants were
%               subjected to the continuous darkness treatment.
%
% Output
%
%       data - structure array containing the fields:
%              data.metNames - List with metabolite names
%              data.abundance - Vector with metabolite abundance
%              data.std - Vector with standard deviation 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1). Verify that absolute values for lipids were retrieved for accessions:
pathMets1 = fullfile('..\ComputeAbsoluteLipids\OutputData\AbsoluteLipids.xlsx');

if ~exist(pathMets1,'file')
    pathFxn = dir('..\ComputeAbsoluteLipids\GetLipidsDarkness.m');
    olDir = cd(pathFxn.folder);
    GetLipidsDarkness % Run the "GetLipidsDarkness" function to obtain the absolute values for lipids
    cd(olDir);
end

% 2). Load data of lipid species whose abudance values are assumed to
% be constant under the conditions evaluated:
pathMets2 = fullfile('InputData','QuantitativeLipidsAccessions.xlsx');
[~, ~, ConstantLipids] = xlsread(pathMets2,'LipidsConstantLevels','A1:D153');
ConstantLipids = RemoveIsNaN(ConstantLipids);

% 2). Load data for chlorophyll abundance:
[~, ~, Chlorophyll] = xlsread(pathMets2,'Chlorophyll','A1:E4');
Chlorophyll = RemoveIsNaN(Chlorophyll);

% 3). Load accessions lipids data:
switch Day
    case 0
        DataSet = char('StandardConditions');
        Range = char('A1:D163');
        ChlID = ['DD', num2str(Day)];
        pathMets = pathMets2;
    case 3
        DataSet = char('AbsoluteLipidsD3');
        Range = char('A1:JZ347');
        ChlID = ['DD', num2str(Day)];
        pathMets = pathMets1;
    case 6
        DataSet = char('AbsoluteLipidsD6');
        Range = char('A1:KB347');
        ChlID = ['DD', num2str(Day)];
        pathMets = pathMets1;
    case {'6h', '21h'}
        DataSet = char('AbsoluteLipidsCol0');
        Range = char('A1:F332');
        ChlID = 'DD0'; % Chlorophyll value is assumed to be equal to the control
        pathMets = pathMets1;
end


% 4). Load accessions data and add metabolites with constant abundance
%     values and clorophyll:
[~, ~, LipidAccessions] = xlsread(pathMets,DataSet,Range);
LipidAccessions = RemoveIsNaN(LipidAccessions);

switch Day
    case 0
        Names = [LipidAccessions(:,1:2); ConstantLipids(2:end,1:2);['Chl',Chlorophyll(strcmp(Chlorophyll(:,1),ChlID),2)]];
    otherwise
        Names = [LipidAccessions(:,1:2);['Chl',Chlorophyll(strcmp(Chlorophyll(:,1),ChlID),2)]];
end
clearvars pathMets idx


% 5). Assign the respective LipidMaps code to the metabolite IDs:

% Load LipidMaps codes:
Names = addLipidMapsCodes(Names);

% 6). Create structure array with the required fields:
DataType = isa(Accession, 'double');
switch DataType
    case 0
        switch Accession
            case {'All'}
                data.metNames  = cellstr(Names(2:end,1));
                Abundances = zeros(size(LipidAccessions,1),size(LipidAccessions,2)-4);
                Abundances(1:size(LipidAccessions,1)-1,:) = str2double(string(LipidAccessions(2:end,5:end)));
                Abundances(end,:) = str2double(string(repmat(Chlorophyll(strcmp(Chlorophyll(:,1),ChlID),3),1,size(Abundances,2))));
                data.Accession = LipidAccessions(1,5:end);
                
            case {'Col0'}
                data.metNames  = cellstr(Names(2:end,1));
                idxAccession = contains(LipidAccessions(1,:), Day);
                Abundances = zeros(size(LipidAccessions,1),1);
                Abundances(1:size(LipidAccessions,1)-1,:) = str2double(string(LipidAccessions(2:end, idxAccession)));
                Abundances(end,:) = str2double(string(Chlorophyll(strcmp(Chlorophyll(:,1),ChlID),3)));
                data.Accession = LipidAccessions(1, idxAccession);
        end
        data.abundance = Abundances;
        SDs = [LipidAccessions(2:end,4);Chlorophyll(strcmp(Chlorophyll(:,1),ChlID),4)];
        data.std = str2double(string(SDs));
    
    case 1
        switch Day
            case 0
                data.metNames  = cellstr(Names(2:end,1));
                Abundances = [LipidAccessions(2:end,2+Accession);ConstantLipids(2:end,3);Chlorophyll(strcmp(Chlorophyll(:,1),ChlID),3)];
                data.abundance = str2double(string(Abundances));
                SDs = [LipidAccessions(2:end,end);ConstantLipids(2:end,4);Chlorophyll(strcmp(Chlorophyll(:,1),ChlID),4)];
                data.std = str2double(string(SDs));
                data.Accession = LipidAccessions(1,3);

            otherwise
                data.metNames  = cellstr(Names(2:end,1));
                Abundances = [LipidAccessions(2:end,4+Accession);Chlorophyll(strcmp(Chlorophyll(:,1),ChlID),3)];
                data.abundance = str2double(string(Abundances));
                SDs = [LipidAccessions(2:end,4);Chlorophyll(strcmp(Chlorophyll(:,1),ChlID),4)];
                data.std = str2double(string(SDs));
                data.Accession = LipidAccessions(1,4+Accession);
        end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%