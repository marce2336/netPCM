function backName = getBackboneName(specName, tailCode)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modified from:
%
% backName = getBackboneName(specName)
%
% Benjamín J. Sánchez (2018-03-24)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize variables:

subspeciesList = {'CDP-DG', 'CL', 'DG', 'DGDG', 'MGDG', 'PA', 'PC', 'PE',...
    'PG', 'PGP', 'PI', 'PS', 'SQDG', 'TG', 'WE'};

%       Lipid ID                Name                       Backbone abbreviation
group1 = {'CL'              'Cardiolipin'                           'CL'
          'CP'              'cutin polymer'                         'CP'
          'Chl'             'Chlorophyll'                           'Chl'
          'DCA'             'dicarboxylic acid'                     'DCA'
          'DG'              'Diacylglycerol'                        'DG'
          'DGDG'            'Digalactosyldiacylglycerol'            'DGDG'
          'F'               'alkyl ferulate'                        'F'
          'FFA'             'free fatty acid'                       'FA'
          'FA'              'acylated fatty acid'                   'FA'
          'FAL'             'fatty aldehyde'                        'FAld'
          'FOH'             'fatty alcohol'                         'FOH'
          'HC'              'Hydrocarbon'                           'HC'
          'Oxo'             'hydroxy fatty acid'                    'OxoFA'
          'HFA'             'hydroxy fatty acid'                    'HFA'
          'MGDG'            'Monogalactosyldiacylglycerol'          'MGDG'
          'OHC'             'Oxygenated hydrocarbon'                'OHC'
          'PC'              'Phosphatidylcholine'                   'PC'
          'PE'              'Phosphatidylethanolamine'              'PE'
          'PG'              'Phosphatidylglycerol'                  'PG'
          'PI'              'Phosphatidylinositol'                  'PI'
          'PS'              'Phosphatidylserine'                    'PS'
          'SQDG'            'Sulfoquinovosyldiacylglycerol'         'SQDG'
          'Campesterol'     'Sterol'                                'Sterol'
          'Sitosterol'      'Sterol'                                'Sterol'
          'Stigmasterol'    'Sterol'                                'Sterol'
          'Brassicasterol'  'Sterol'                                'Sterol'
          'TG'              'Triacylglycerol'                       'TG'
          'WE'              'Wax esther'                            'Wax'
          'Cer'             'Ceramide'                              'Cer'
          'GlcCer'          'glycosylceramide'                      'GCer'
          'Gal-GlcA-IPC'    'GIPC'                                  'Hex_GlcA_IPC'
          'Man-GlcA-IPC'    'GIPC'                                  'Hex_GlcA_IPC'
          'Glc-GlcA-IPC'    'GIPC'                                  'Hex_GlcA_IPC'
          'GlcA-IPC'        'GIPC'                                  'GlcA_IPC'
          'Ara-GlcA-IPC'    'GIPC'                                  'Pen-GlcA_IPC'
          'CDP-DG'          'CDP-diacylglycerol'                    'CDP-DG'
          'PA'              'Phosphatidic acid'                     'PA'
          'PGP'             'Diacylglycerophosphoglycerophosphate'  'PGP'
          'MGMG'            'Monogalactosylmonoacylglycerol'        'MGMG'
          'LPC'             'Lyso-phosphatidylcholine'              'LPC'
          'LPA'             'Lyso-phosphatidic acid'                'LPA'};


backName = '';

% Find the group to which the backbone belongs:
Flag =  sum(strcmp(subspeciesList, specName));

switch Flag
    case 0
        backName = RetrieveBackbone(group1, specName);
        backName = [backName '-backbone'];
        
    case 1
        if exist('tailCode', 'var')
            backName = RetrieveBackbone(group1, specName);
            tailCode = strsplit(tailCode, ':');
            backName = [backName,tailCode{1},'-backbone'];
        else
            backName = RetrieveBackbone(group1, specName);
            backName = [backName '-backbone'];
        end
end

end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function backName = RetrieveBackbone(group1, specName)

for i = 1:length(group1)
    if startsWith(specName,group1{i,1})
        backName = group1{i,3};
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%