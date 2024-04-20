function TailsModel = ObtainModelTails(backCode, MatchedMets)

% Initialize variables:
Cols = length(find(char(MatchedMets(:,2)) == ':'));
TailsModel = cell(size(MatchedMets,1),Cols/size(MatchedMets,1));

% Obtain acyl chains from matched lipids:
for i = 1:size(MatchedMets,1)
    Tail = MatchedMets{i,2};
    
    switch backCode
        case {'PS', 'PE', 'PC', 'PI', 'SQDG', 'PG', 'DG', 'MGDG', 'DGDG', 'TG',...
                'CDP-DG', 'LPA', 'LPC', 'MGMG', 'PA', 'PGP', 'WE'}
            TailsModel(i,:) = strsplit(Tail,{'/','_'});
            
        case {'CL'}
            TailsModel(i,:) = strsplit(Tail,{'/', '_', ',3''-['});
    end
end

% Eliminate structural information:
switch backCode
    case {'PS', 'PE', 'PC', 'PI', 'SQDG', 'PG', 'DG', 'MGDG', 'DGDG','CDP-DG', ...
            'PA', 'PGP', 'WE'}
        TailsModel(:,1) = extractBetween(TailsModel(:,1), 2, 5);
        TailsModel(:,2) = extractBetween(TailsModel(:,2), 1, 4);
            
    case {'TG'}
        TailsModel(:,1) = extractBetween(TailsModel(:,1), 2, 5);
        TailsModel(:,2) = extractBetween(TailsModel(:,2), 1, 4);
        TailsModel(:,3) = extractBetween(TailsModel(:,3), 1, 4);
        
    case {'LPA', 'LPC'}
        % Examples: '[GP1005] LPA(16:0/0:0)'
        %           '[GP0105] LPC(16:0/0:0)'
        TailsModel(:,1) = extractBetween(TailsModel(:,1), 2, 5);
        TailsModel(:,2) = extractBetween(TailsModel(:,2), 1, 3);
        
    case {'MGMG'}
        % Example: '[GL0401] MGMG(0:0/16:2(7Z,10Z))'
        TailsModel(:,1) = extractBetween(TailsModel(:,1), 2, 4);
        TailsModel(:,2) = extractBetween(TailsModel(:,2), 1, 4);
        
    case {'CL'}
        TailsModel(:,1) = strrep(TailsModel(:,1), '1''-[', '');
        TailsModel(:,1) = extractBetween(TailsModel(:,1), 2, 5);
        TailsModel(:,2) = extractBetween(TailsModel(:,2), 1, 4);
        TailsModel(:,3) = extractBetween(TailsModel(:,3), 1, 4);
        TailsModel(:,4) = extractBetween(TailsModel(:,4), 1, 4);    
end

end