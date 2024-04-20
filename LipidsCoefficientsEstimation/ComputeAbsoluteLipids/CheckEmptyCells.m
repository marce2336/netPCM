function missingAbsolute = CheckEmptyCells(missingAbsolute, abundanceValues, missingNames, absoluteNames)

% Find species that were not detected:
isEmpty = sum((isnan(missingAbsolute)),2);
uniqueNames = missingNames((isEmpty ~= 0), 1);
uniqueNames = unique(uniqueNames);

% Find the respective class of missing lipids:
for i = 1:numel(uniqueNames)
    lipid_i = uniqueNames{i};

    switch lipid_i
        case {'PC','PE','PG','PI','PS','PA','LPA','LPC','PGP'} % Phospholipids
            species = {'PC','PE','PG','PI','PS','PA','LPA','LPC','PGP'};
            absoluteClass_i = GetValuesForClass(absoluteNames, abundanceValues, species);

        case {'MG','DG','TG','MGMG','MGDG','DGDG','CL','SQDG','CDP-DG','WE'} % Neutral lipids
            species = {'MG','DG','TG','MGMG','MGDG','DGDG','CL','SQDG','CDP-DG','WE'};
            absoluteClass_i = GetValuesForClass(absoluteNames, abundanceValues, species);

        case {'FFA','FA','FAL','FOH','HC','OHC','Oxo','DCA','HFA'} % Fatty acyls
            species = {'FFA','FA','FAL','FOH','HC','OHC','Oxo','DCA','HFA'};
            absoluteClass_i = GetValuesForClass(absoluteNames, abundanceValues, species);

        case {'SPB','Cer','GlcCer','GlcA-IPC','Ara-GlcA-IPC','Gal-GlcA-IPC','Glc-GlcA-IPC','Man-GlcA-IPC','PI-Cer'} % Sphingolipids
            species = {'SPB','Cer','GlcCer','GlcA-IPC','Ara-GlcA-IPC','Gal-GlcA-IPC','Glc-GlcA-IPC','Man-GlcA-IPC','PI-Cer'};
            absoluteClass_i = GetValuesForClass(absoluteNames, abundanceValues, species);

        case {'Campesterol','Sitosterol','Stigmasterol','ST','Brassicasterol'} % Sterols
            species = {'Campesterol','Sitosterol','Stigmasterol','ST','Brassicasterol'};
            absoluteClass_i = GetValuesForClass(absoluteNames, abundanceValues, species);

    end

    % Find indexes for missing data and fill according to respective class:
    idxClass = ismember(missingNames(:,1), lipid_i);
    missingAbsolute(idxClass,:) = absoluteClass_i;

end

end
