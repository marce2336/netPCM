function newNames = addLipidMapsCodes(names)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The 'addLipidMapsCodes' function retrieves trhe respective Lipid Maps
% codes for a list of lipid species provided. Each lipid species must be
% identified by the lipid abrreviation, followed by the sum composition, as
% indicated in the example next:
%
% 'Abbreviation'   'Total carbon atoms'   ':'   'total unsaturations'
%  TG 52:3
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load LipidMaps codes:
pathLMaps = fullfile('InputData','LipidMapsClassification.xlsx');
[~, ~, LipidMaps] = xlsread(pathLMaps,'LM_Classification','A1:E45');
LipidMaps = RemoveIsNaN(LipidMaps);

% Add LipidMaps codes to metabolite names:
for i = 2:size(names(:,1),1)
    LipidID = names(i,1);
    splitID = strsplit(LipidID{1,1},' ');
    
    switch splitID{1,1}
        case {'Cer'}
            CerTail = strsplit(splitID{1,2}, '/');
            
            if contains(CerTail{1,1}, 'd18:1') || contains(CerTail{1,1}, 't18:1')
                LMidx = contains(LipidMaps(:,1), splitID{1,1}).*contains(LipidMaps(:,1),CerTail{1,1}) == 1;
                
            elseif contains(CerTail{1,1}, 'OH')
                if length(strfind(CerTail{1,1}, 'OH')) == 2
                    TailID = ['d', extractBefore(CerTail{1,1},';')];
                    
                elseif length(strfind(CerTail{1,1}, 'OH')) == 3
                    TailID = ['t', extractBefore(CerTail{1,1},':')];
                end
                
                LMidx = contains(LipidMaps(:,1), splitID{1,1}).*contains(LipidMaps(:,1),TailID) == 1;
            end
            
        otherwise
            LMidx = strcmp(LipidMaps(:,1), splitID{1,1});
    end
    
    findCode = find(~cellfun(@isempty, LipidMaps(LMidx,:)));
    LMcode = extractAfter(LipidMaps(LMidx,size(findCode,2)), '[');
    
    if size(names, 2) > 1
        LipidStructure = names{i,2};
    else
        LipidStructure = '';
    end
    
    if isempty(LipidStructure)
        names{i,1} = ['[', LMcode{1,1}, ' ', LipidID{1,1}];
    else
        names{i,1} = ['[', LMcode{1,1}, ' ', LipidID{1,1}, ' ', LipidStructure];
    end
end

newNames = names;

end

