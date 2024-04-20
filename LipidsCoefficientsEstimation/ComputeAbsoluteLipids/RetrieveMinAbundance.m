function minAbundance = RetrieveMinAbundance(LipidNames, abundanceValues, flag)

% Create array with the names of lipid classes included in the list provided:
uniqueNames  = unique(LipidNames(:,1));

% Search for missing classes in the unique names list and act accordingly:
uniqueNames = VerifyMetsClasses(uniqueNames);

switch flag
    case 1 % Retrieve minimum value of abundance for dataset of treated samples:
        minAbundance = zeros(numel(uniqueNames),1);

        for i = 1:numel(uniqueNames)
            lipidClass_i = uniqueNames{i};
            classIdx = strcmp(LipidNames(:,1), lipidClass_i);
            minLipidClass_i = min(min((abundanceValues(classIdx, :))));

            % If there is not a single species of the category, find the
            % minimum value existing in the data set:
            if isempty(minLipidClass_i)
                minLipidClass_i = min(min(str2double(string(abundanceValues(:, :)))));
            end

            minAbundance(i) = minLipidClass_i;
        end

    case 2 % Retrieve minimum value of abundance for dataset of control samples:
        minAbundance = zeros(numel(uniqueNames),2);

        for i = 1:numel(uniqueNames)
            lipidClass_i = uniqueNames{i};
            classIdx = strcmp(LipidNames(:,1), lipidClass_i);
            minLipidClass_i = min(min(str2double(string(abundanceValues(classIdx, 1)))));
            meanSDclass_i = min(min(str2double(string(abundanceValues(classIdx, 2)))));

            % If there is not a single species of the category, find the
            % minimum value existing in the data set:
            if isempty(minLipidClass_i)
                minLipidClass_i = min(abundanceValues(:,1));
                meanSDclass_i   = min(abundanceValues(:,2));
            end

            minAbundance(i,1) = minLipidClass_i;
            minAbundance(i,2) = meanSDclass_i; 
        end
end

minAbundance = cellstr(string(minAbundance));
minAbundance = [uniqueNames, minAbundance];

end