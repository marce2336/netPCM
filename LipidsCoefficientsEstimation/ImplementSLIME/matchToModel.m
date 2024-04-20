function [pos, backbone, chains] = matchToModel(model,metName, flagIdxs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modified from:
%
% pos = matchToModel(model,metName)
%
% Benjamín J. Sánchez. Last update (2018-03-25)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Translating codes to names in model:
%codes = {'TG'               'triacylglycerol'                         '[GL0301]'
%         'DG'               'diacylglycerol'                          '[GL0201]'
%         'MG'               'monoacylglycerol'                        '[GL0101]'
%         'PS'               'phosphatidyl-L-serine'                   '[GP0301]'
%         'PE'               'phosphatidylethanolamine'                '[GP0201]'
%         'PC'               'phosphatidylcholine'                     '[GP0101]'
%         'PI'               '1-phosphatidyl-1D-myo-inositol'          '[GP0601]'
%         'PG'               'phosphatidylglycerol'                    '[GP0401]'
%         'CL'               'cardiolipin'                             '[GP1201]'
%         'MGDG'             'monogalactosyldiacylglycerol'            '[GL0501]'
%         'DGDG'             'digalactosyldiacylglycerol'              '[GL0501]'
%         'SQDG'             'sulfoquinovosyldiacylglycerol'           '[GL0501]'
%         'CDP-DG'           'CDP-diacylglycerol'                      '[GP1301]'
%         'LPA'              'monoacylglycerophosphate'                '[GP1005]'
%         'LPC'              'Monoacylglycerophosphocholines'          '[GP0105]'
%         'MGMG'             'monogalactosylmonoacylglycerol'          '[GL0401]'
%         'PA'               'diacylglycerophosphate'                  '[GP1001]'
%         'PGP'              'diacylglycerophosphoglycerophosphate'    '[GP0501]'
%         'WE'               'wax monoesters'                          '[FA0701]'
%         'FFA'              'free fatty acyl'                         '[FA0101]'
%         'FA'               'fatty acyl'                              '[FA0101]'
%         'FAL'              'fatty aldehyde'                          '[FA06]'
%         'FOH'              'fatty alcohol'                           '[FA05]'
%         'OHC'              'Oxygenated hydrocarbon'                  '[FA12]'
%         'HC'               'Hydrocarbons'                            '[FA11]'
%         'Oxo'              'Oxo fatty acid'                          '[FA0106]'
%         'DCA'              'Dicarboxylic acid'                       '[FA0117]'
%         'HFA'              'Hydroxy fatty acid'                      '[FA0105]'
%         'Chl'              'Chlorophyll'                             '[PR010401]'
%         'Campesterol'      'Campesterol'                             '[ST0103]'
%         'Sitosterol'       'Sitosterol'                              '[ST0104]'
%         'Stigmasterol'     'Stigmasterol'                            '[ST0104]'
%         'Brassicasterol'   'Brassicasterol'                          '[ST0103]'};  
%----------------------------------------------------------------------------------

if ~exist('flagIdxs', 'var')
    flagIdxs = false;
end

pos = false(size(model.mets));

%Split metabolite into backbone and tails:
parts    = strsplit(metName,' ');
if length(parts) == 2
    LMcode = parts{1};
    backCode = parts{2};
elseif length(parts) == 3
    LMcode = parts{1};
    backCode = parts{2};
    tailCode = parts{3};
elseif length(parts) == 4
    LMcode = parts{1};
    backCode = parts{2};
    tailCode = parts{3};
    tailStructure = parts{4};
end

if contains(LMcode,'[PR010401]') || contains(LMcode,'[ST0103]') ||...
        contains(LMcode,'[ST0104]')
    
    %Find mets with matching LipidMaps code:
    findMet = contains(model.metNames,LMcode).*contains(model.metNames,backCode) == 1;
    
    if sum(findMet) > 0
        pos(findMet) = true;
        specName = backCode;
        backbone = getBackboneName(specName);
        
        switch LMcode
            case {'[PR010401]'}
                chains = 'Phytyl';
                
            case {'[ST0103]', '[ST0104]'}
                chains = backCode;
        end
        
    else
        backbone = '';
        chains = '';
    end
    
else    
              
    switch backCode
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Lipid species with one (1) or more fatty acyl chains
        case {'PS', 'PE', 'PC', 'PI', 'SQDG', 'PG','DG', 'MGDG', 'DGDG', 'TG',...
                'CDP-DG', 'CL', 'LPA', 'LPC', 'MGMG', 'PA', 'PGP', 'WE'}
            
            %If only sum composition is provided, find mass isomers and act accordingly:
            if length(parts) == 3 
                Hits = contains(model.metNames,LMcode).*contains(model.metNames,backCode) == 1;
                [MatchedMets, findMet] = MatchSumComposition(Hits, model, tailCode, backCode);
            
            %If structure data is provided, find tails and act accordingly:
            elseif length(parts) == 4
                Hits = contains(model.metNames,LMcode).*contains(model.metNames,backCode) == 1;
                [MatchedMetCodes, findMetIdx] = MatchSumComposition(Hits, model, tailCode, backCode);
                [findMet,MatchedMets] = CompareTailStructure(tailStructure, MatchedMetCodes, findMetIdx, backCode);
            end
            
            if sum(sum(~cellfun(@isempty, MatchedMets))) > 0
                % Obtain tails from matched lipids in model:
                TailsModel = ObtainModelTails(backCode, MatchedMets);
                
                % Obtain backbone name:
                specName = backCode;
                backbone = getBackboneName(specName, tailCode);
                
                % Sort chains to eliminate isomers with identical acyl chains:
                %Function downloaded from Mathworks: Jim Hokanson (2021). Unique Rows for a cell array
                %(https://www.mathworks.com/matlabcentral/fileexchange/25917-unique-rows-for-a-cell-array), 
                %MATLAB Central File Exchange. Retrieved August 18, 2021.
                [chains, IsLipid] = uniqueRowsCA(TailsModel);
                
                % Obtain indexes of mass isomers in model:
                Idxs = findMet(IsLipid);
                
                switch flagIdxs
                    case false
                        pos = unique(Idxs);     
                    case true
                        pos = findMet;
                end


            else
                backbone = '';
                chains = '';
            end
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Identification of features for fatty acyls:   
        case {'FFA', 'FA', 'FAL', 'HC', 'FOH', 'OHC', 'Oxo', 'DCA', 'HFA'}
            % Split tail:
            SplitTail = strsplit(tailCode,{':',';'});
            
            %Find equivalent compounds in model:
            Hits = contains(model.metNames,LMcode);
            Idxs = find(Hits);
            MatchedMets = model.metNames(Hits);

            if sum(sum(~cellfun(@isempty, MatchedMets))) > 0
                % Eliminate Fatty acyl-CoAs from list:
                IsCoA = contains(MatchedMets(:,1),cellstr('[FA0705]'));
                MatchedMets(IsCoA == 1, :) = '';
                Idxs(IsCoA == 1) = '';
                MatchedMets = split(MatchedMets,' ');

                %Find matching metabolites in model:

                % CASE 1: the fatty acyl contains functional groups and has
                %  full information at the structure level.
                if size(SplitTail,2) >= 3
                    IsLipid = find(strcmp(tailCode,MatchedMets(:,3)));

                % CASE 2: The fatty acyl has no functional groups and has
                % full information at the structure level.
                elseif size(SplitTail,2) == 2 && strlength(SplitTail{1,2}) > 2
                    IsLipid = find(strcmp(tailCode,MatchedMets(:,3)));

                % CASE 3: The fatty acyl has no functional groups and there is
                % information available only at the species level.    
                elseif size(SplitTail,2) == 2 && strlength(SplitTail{1,2}) < 2
                    IsLipid = contains(string(MatchedMets(:,3)),tailCode);
                end

                if sum(IsLipid) > 0
                    specName = backCode;
                    backbone = getBackboneName(specName);
                    pos = Idxs(IsLipid);

                 switch LMcode
                     case {'[FA0101]', '[FA11]'}
                         chains = tailCode;
                     case {'[FA06]'}
                         chains = ['FAL-',tailCode];
                     case {'[FA05]'}
                         chains = ['FOH-',extractBefore(tailCode,';')];
                     case {'[FA12]'}
                         chains = ['C',SplitTail{1,1},'-',extractBefore(SplitTail{1,3},'O'),'k'];
                     case {'[FA0106]'}
                         if strlength(SplitTail{1,2}) > 2
                             DB = extractBefore(SplitTail{1,2},'(');
                         else
                             DB = SplitTail{1,2};
                         end
                         chainDB = str2double(string(DB))-1;
                         chains = ['C',SplitTail{1,1},':',num2str(chainDB),'-w',SplitTail{1,3}];
                     case {'[FA0117]'}
                         if strlength(SplitTail{1,2}) > 2
                             DB = extractBefore(SplitTail{1,2},'(');
                         else
                             DB = SplitTail{1,2};
                         end
                         chains = ['C',SplitTail{1,1},':',DB,'-DCA'];
                     case {'[FA0105]'}
                         if length(strfind(tailCode,'OH')) == 2
                             chains = ['C',extractBefore(tailCode,';'),'-',...
                                 extractBefore(SplitTail{1,3},'OH'),'-',extractBefore(SplitTail{1,4},'OH'),'diOH']; %C18-triOH
                         elseif length(strfind(tailCode,'OH')) == 3
                             chains = ['C',SplitTail{1,1},'-triOH'];
                         end
                 end

                else
                    backbone = '';
                    chains = '';
                end

            else
                backbone = '';
                chains = '';
            end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Identification of features for sphingolipids:
        case {'Cer', 'GlcCer', 'GlcA-IPC', 'Gal-GlcA-IPC', 'Man-GlcA-IPC', 'Glc-GlcA-IPC', 'Ara-GlcA-IPC'}
            %Find similar compounds in model:
            Hits = contains(model.metNames,LMcode);
            Idxs = find(Hits);
            MatchedMets = model.metNames(Hits);

            if sum(sum(~cellfun(@isempty, MatchedMets))) > 0
                MatchedMets = split(MatchedMets,{' ','/'});
                Tails = strsplit(tailCode, '/');
                LCB = strsplit(Tails{1,1}, ';'); % Get information from long chain base
                ShortTail = strsplit(Tails{1,2}, ';');

                % CASE 1: full information at the structure level is provided for the LCB and the short chain does not contain OH-groups.
                if size(LCB,2) > 1 && size(ShortTail,2) <= 1
                    sphingoVariant = 1;
                    [Idxs, chains] = MatchSphingoVariant(sphingoVariant, MatchedMets, LCB, ShortTail, Idxs, backCode);
                    
                % CASE 2: full information at the structure level is provided for the LCB and the short chain contain OH-groups.
                elseif size(LCB,2) > 1 && size(ShortTail,2) > 1
                    sphingoVariant = 2;
                    [Idxs, chains] = MatchSphingoVariant(sphingoVariant, MatchedMets, LCB, ShortTail, Idxs, backCode);

                % CASE 3: information at the species level is provided for the LCB and the short chain contain OH-groups.
                elseif size(LCB,2) <= 1 && size(ShortTail,2) > 1
                    sphingoVariant = 3;
                    [Idxs, chains] = MatchSphingoVariant(sphingoVariant, MatchedMets, LCB, ShortTail, Idxs, backCode);
                end

                if sum(Idxs) > 0
                    specName = backCode;
                    backbone = getBackboneName(specName);
                    pos = Idxs;

                else
                    backbone = '';
                    chains = '';
                end

            else
                backbone = '';
                chains = '';
            end
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
        % Identification of features for glycerolipids containing fatty acyls with functional groups:
        case {'MG'}
            %Find equivalent compounds in model:
            Hits = contains(model.metNames,LMcode) == 1;
            Idxs = find(Hits);
            MatchedMets = model.metNames(Hits);
            
            if sum(sum(~cellfun(@isempty, MatchedMets))) > 0
                MatchedMets = split(MatchedMets,backCode);
                SplitTail = strsplit(tailStructure,'/');
                SplitTail{1,1} = strrep(SplitTail{1,1},'(', '');
                SplitTail{1,size(SplitTail,2)} = strrep(SplitTail{1,size(SplitTail,2)},')', '');
                IsLipid = contains(MatchedMets(:,2),(SplitTail{1,1})).*contains(MatchedMets(:,2),SplitTail{1,2})...
                    .*contains(MatchedMets(:,2),(SplitTail{1,3})) == 1;

                if sum(IsLipid) > 0
                    specName = backCode;
                    backbone = getBackboneName(specName);
                    pos = Idxs(IsLipid);
                    SplitTail = SplitTail(~contains(SplitTail, '0:0') == 1);
                    chains = cell(1, size(SplitTail,2));

                    for k = 1:size(SplitTail,2)
                        Tail = SplitTail(1,k);
                        TailComponents = split(Tail,';');

                        if sum(contains(TailComponents,'OH')) == 2
                            chainID = ['C',TailComponents{1,1},'-',extractBefore(TailComponents{2,1},'OH'),...
                                '-',extractBefore(TailComponents{3,1},'OH'),'diOH'];
                            chains{k} = chainID;
                        elseif sum(contains(TailComponents,'OH')) == 3
                            chainID = ['C',extractBefore(TailComponents{1,1},':'),'-triOH'];
                            chains{k} = chainID;
                        elseif sum(contains(TailComponents,'O2')) >= 1
                            chainID = ['C',tailCode,'-DCA'];
                            chains{k} = chainID;
                        end
                    end
                    
                else
                    backbone = '';
                    chains = '';
                end
                
            else
                backbone = '';
                chains = '';
            end
    end
end
end
