function GetLipidsDarkness

%% Estimation of absolute lipids content for Arabidopsis accessions:
%  This function allows the estimation of the absolute lipids content for  
%  the Arabidopsis accessions, using quantitative data from plants grown 
%  under standard conditions retrieved from Liu et al. (2020) Mol.
%  Plant. 13:1523–1532 (https://doi.org/10.1016/j.molp.2020.07.016), and
%  Li-Beisson et al. (2010) Arabidopsis Book 8:e0133 (doi:
%  10.1199/tab.0133); and the relative abundance data obtained for all the
%  accessions. To this end, first the ratios of Dark- vs. Light-grown 
%  plants are calculated using the relative abundance data by implementing
%  the "CalculateRatios" function. Then, the metabolites content is 
%  adjusted for each accession by multiplying the quantitative data for 
%  plants grown under standard conditions by its respective abundance
%  ratios.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1). Load quantitative data of lipids under standard conditions:
%      The quantitative data includes the following lipid classes: PC, PE,
%      PG, PI, PS, DG, TG, MGDG, DGDG.

pathMets = fullfile('InputData','LipidsCompiled.xlsx');
[~, ~, QLipidsControl] = xlsread(pathMets,'MergedData_g_g-1DW','A1:D267');
QLipidsControl = RemoveIsNaN(QLipidsControl);

%% 2). Load ratios obtained from relative abundance data of lipids under extended darkness:
%  NOTE: only lipids present in the LipidModule are included in the list

pathMets = fullfile('..\CalculateAbsoluteLipids\OutputData\Ratios.xlsx');

if ~exist(pathMets,'file')
    CalculateRatios % Run the "CalculateRatios" function to obtain the lipid ratios:
end

pathMets = fullfile('OutputData','Ratios.xlsx');
[~, ~, LipidRatiosD3] = xlsread(pathMets,'RatiosD3','A1:HX285');
LipidRatiosD3 = RemoveIsNaN(LipidRatiosD3, 1);

[~, ~, LipidRatiosD6] = xlsread(pathMets,'RatiosD6','A1:HZ285');
LipidRatiosD6 = RemoveIsNaN(LipidRatiosD6, 1);

[~, ~, LipidRatiosCol0] = xlsread(pathMets,'RatiosCol0','A1:FF3');
LipidRatiosCol0 = RemoveIsNaN(LipidRatiosCol0, 2);

%% 3). Match lipid species under extended darkness (D3) and standard conditions:
QLipidsD3 = MatchLipidSpecies(QLipidsControl, LipidRatiosD3);


%% 4). Match lipid species under extended darkness (D6) and standard conditions:
QLipidsD6 = MatchLipidSpecies(QLipidsControl, LipidRatiosD6);


%% 5). Match lipid species of Col-0 under extended darkness (6h and 21h) and standard conditions:
QLipidsCol0 = MatchLipidSpecies(QLipidsControl, LipidRatiosCol0);

%% 6). Save output variables:
pathMets = fullfile('OutputData','AbsoluteLipids.xlsx');
xlswrite(pathMets, QLipidsD3, 'AbsoluteLipidsD3' );
xlswrite(pathMets, QLipidsD6, 'AbsoluteLipidsD6');
xlswrite(pathMets, QLipidsCol0, 'AbsoluteLipidsCol0');

end