function missingAbsoluteNew = LipidsUnderDetection(missingMets, absoluteLipids, flag, lipidRatios)

absoluteNames    = SplitNames(absoluteLipids(2:end,1));
missingNames     = SplitNames(missingMets);
uniqueMissing     = unique(missingNames(:,1));

switch flag
    case 1 % Lipids are detected in the control but not in treated samples:
        % Calculate minimum value for lipid classes in treated samples:
        missingAbsolute  = zeros(numel(missingMets),size(absoluteLipids,2)-4);
        abundanceValues = absoluteLipids(2:end,5:end);
        isEmpty = cellfun(@isempty, abundanceValues);
        abundanceValues(isEmpty) = {' '};
        abundanceValues = str2double(string(abundanceValues));
        minAbundance = RetrieveMinAbundance(absoluteNames, abundanceValues, 1);

        % Fill missing values in treated samples with 1/5 of the minimum
        % positive value detected for the respective lipid class:
        minAbundanceValues = str2double(string(minAbundance(:,2)));
        minAbundanceValues = minAbundanceValues./5;

        for i = 1:numel(uniqueMissing)
            missing_i = uniqueMissing{i};
            idxMissingList = strcmp(missingNames(:,1), missing_i);
            idxMinValue_i = strcmp(minAbundance(:,1), missing_i);
            missingAbsolute(idxMissingList,:) = minAbundanceValues(idxMinValue_i);
        end

        % Search for blank spaces corresponding to lipid classes for which
        % none species was detected and act accordingly:
        missingAbsoluteNew = CheckEmptyCells(missingAbsolute, abundanceValues, missingNames, absoluteNames);

    case 2 % Lipids are detected in treated samples but not in the control:
        missingAbsolute  = zeros(numel(missingMets),size(absoluteLipids,2)-2);

        % Calculate minimum value for lipid classes in the control:
        abundanceValues = str2double(string(absoluteLipids(2:end,3:4)));
        minAbundance = RetrieveMinAbundance(absoluteNames, abundanceValues, 2);

        % Fill missing values in treated samples with 1/5 of the minimum
        % positive value detected for the respective lipid class:
        minAbundanceValues = str2double(string(minAbundance(:,2)));
        minAbundanceValues = minAbundanceValues./5;
        minSDvalues = str2double(string(minAbundance(:,3)));
        minSDvalues = minSDvalues./8;

        for i = 1:numel(uniqueMissing)
            missing_i = uniqueMissing{i};
            idxMissingList = strcmp(missingNames(:,1), missing_i);
            idxMinValue_i = strcmp(minAbundance(:,1), missing_i);
            missingAbsolute(idxMissingList,1) = minAbundanceValues(idxMinValue_i);
            missingAbsolute(idxMissingList,2) = minSDvalues(idxMinValue_i);
        end

        % Adjust abundance taking into account fold change:
        missingAbsoluteNew = AdjustFoldChange(missingMets, missingAbsolute, lipidRatios, absoluteLipids);
end

end