function LIPIDcoefficients = calculateLipidsCoeff(Accession,Day)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% modified from function: data = readEjsingData(i), created by Benjamín J.
% Sánchez (Last update: 2018-03-20).
%
% Usage
%
%       LIPIDcoefficients = CalculateLipidsCoeff(Accession,Day)
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
%              LIPIDcoefficients.Accessions - List of accession names
%              LIPIDcoefficients.ChainNames - List of names of lipid chains
%              LIPIDcoefficients.ChainAbundance - Vector of lipid chains abundance
%              LIPIDcoefficients.BackboneNames - Vector of lipid backbone names
%              LIPIDcoefficients.BackboneAbundance - Vector of lipid backbone abundance
%              LIPIDcoefficients.TotalLipid_umol - Total lipids content umol g-1DW
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 1). Load accessions lipids data:

% Load lipids abundance data:
data = readLipidsData(Accession,Day);

% Load metabolic model:
modelPath = fullfile('InputData', 'LipidModule.mat');
model = readCbModel(modelPath);


%% 2). Get list of metabolites participating in lipid pseudoreactions:

[ChainMets, BackboneMets] = obtainPseudometabolites(model);

%% 3). Initialize variables:

% Create array struct with the fields: 
LIPIDcoefficients.Accessions = cell(size(data.Accession)); % Accessions name
LIPIDcoefficients.ChainNames = cell(size(ChainMets)); % Chain names
LIPIDcoefficients.ChainAbundance = zeros(size(ChainMets,1),size(data.Accession,2)); % Chains abundance
LIPIDcoefficients.BackboneNames = cell(size(BackboneMets)); % Backbone names
LIPIDcoefficients.BackboneAbundance = zeros(size(BackboneMets,1),size(data.Accession,2)); % Backbones abundance
LIPIDcoefficients.TotalLipid_umol = zeros(1,size(data.Accession,2)); % Total lipids content umol g-1DW

%% 4). Calculate the SLIMEr coefficients:

DataType = isa(Accession, 'double');
switch DataType
    
    case 1 % Calculate SLIMEr coefficients for the specific accession:
        SLIMEcoeff = GetSLIMErCoeff(data,model,1);
        [ChainCoeff, BackboneCoeff] = getAccessionCoefficients(ChainMets, BackboneMets, SLIMEcoeff);
        AccessionName = num2str(Accession);
        
           % Add data to fields:
        if Day == 0
            LIPIDcoefficients.Accessions = {'Control'};
         else
             LIPIDcoefficients.Accessions = SLIMEcoeff.lipidData.Accession;
        end
        
        LIPIDcoefficients.ChainNames = ChainMets;
        LIPIDcoefficients.ChainAbundance = ChainCoeff;
        LIPIDcoefficients.BackboneNames = BackboneMets;
        LIPIDcoefficients.BackboneAbundance = BackboneCoeff;
        LIPIDcoefficients.TotalLipid_umol = SLIMEcoeff.lipidData.TotalLipid_umol;
        
    case 0 % Calculate SLIMEr coefficients for all accessions:
        
        LIPIDcoefficients.ChainNames = ChainMets;
        LIPIDcoefficients.BackboneNames = BackboneMets;
        AccessionName = Accession;
 
        for i = 1:size(data.Accession,2)
            dataAccession.Accession = data.Accession(1,i);
            dataAccession.metNames = data.metNames;
            dataAccession.abundance = data.abundance(:,i);
            dataAccession.std = data.std;
            SLIMEcoeff = GetSLIMErCoeff(dataAccession,model,1);
            [ChainCoeff, BackboneCoeff] = getAccessionCoefficients(ChainMets, BackboneMets, SLIMEcoeff);
            LIPIDcoefficients.Accessions(1,i) = data.Accession(1,i);
            LIPIDcoefficients.ChainAbundance(:,i) = ChainCoeff;
            LIPIDcoefficients.BackboneAbundance(:,i) = BackboneCoeff;
            LIPIDcoefficients.TotalLipid_umol(1,i) = SLIMEcoeff.lipidData.TotalLipid_umol;
        end         
end

% Eliminate columns with zero entries:
LIPIDcoefficients.ChainAbundance( :, all( ~any( LIPIDcoefficients.ChainAbundance ), 1 ) ) = [];
LIPIDcoefficients.BackboneAbundance( :, all( ~any( LIPIDcoefficients.BackboneAbundance ), 1 ) ) = [];
LIPIDcoefficients.TotalLipid_umol(LIPIDcoefficients.TotalLipid_umol == 0) = [];
LIPIDcoefficients.Accessions = LIPIDcoefficients.Accessions(~cellfun(@isempty, LIPIDcoefficients.Accessions));

% Calculate backbones abundance for lipid pools which contain subspecies:
LIPIDcoefficients = getPoolAbundance(LIPIDcoefficients);

% Save accessions coefficients in output folder:
fileName = ['LIPIDcoefficients','_Day',num2str(Day),'_','Accession','-',AccessionName];
save(fileName, 'LIPIDcoefficients')
pathMets = fullfile('OutputData');
movefile([fileName,'.mat'],pathMets)

end
