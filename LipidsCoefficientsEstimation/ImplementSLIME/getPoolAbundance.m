function LIPIDcoefficients = getPoolAbundance(LIPIDcoefficients)

subspeciesList  = {'CDP-DG', 'CL', 'DG', 'DGDG', 'MGDG', 'PA', 'PC', 'PE',...
    'PG', 'PGP', 'PI', 'PS', 'SQDG', 'TG', 'WE'};

backbones   = LIPIDcoefficients.BackboneNames;
backbones   = extractBefore(backbones, 'backbone');
        
for i = 1:size(LIPIDcoefficients.BackboneNames,1)
    backName = LIPIDcoefficients.BackboneNames{i};
    backName = extractBefore(backName, '-backbone');
    isPool = regexp(backName, '[A-Z]', 'match');
    
    if sum(~cellfun(@isempty, isPool)) > 0 && sum(strcmp(subspeciesList, backName)) > 0
        % Find regular expression:
        pattern     = [backName '+[0-9]'];
        matchSubspecie = regexp(backbones, pattern, 'match');
        subSpecies  = (~cellfun(@isempty, matchSubspecie)).*(strlength(backName)+3 == strlength(backbones)).*(contains(backbones, backName)) == 1;
        
        % Calculate pool abundance:
        abundance = sum(LIPIDcoefficients.BackboneAbundance(subSpecies,:),1);
        LIPIDcoefficients.BackboneAbundance(i,:) = abundance;
    end
end

end
