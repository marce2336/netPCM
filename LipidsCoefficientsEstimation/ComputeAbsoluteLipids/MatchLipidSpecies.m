function absoluteLipids = MatchLipidSpecies(QLipidsControl, lipidRatios)

% ###
% ##### 1). Match information to find common lipid species among treated plants and standard conditions:
% ###
absoluteLipids = [QLipidsControl, cell(size(QLipidsControl,1),size(lipidRatios,1)-1)];
absoluteLipids(1,3) = {'g g–1 DW Control'};
absoluteLipids(1,5:end) = lipidRatios(2:end,1);
%includedInRatios = cell(size(lipidRatios,1),1);
%includedInRatios{1} = string(strcat('Day_', string(day)));
count = 1;
includedControl = zeros(size(QLipidsControl,1),1);

for i = 2:size(QLipidsControl,1)
    lipid_i = QLipidsControl(i,1);
    
    % Get abundance for sample under standard conditions
    controlAbsoluteQ = str2double(string(QLipidsControl(i,3)));
    [~,col_i] = find(strcmp((lipidRatios(1,1:end)),lipid_i));
    
    if ~isempty(col_i)
        count = count + 1;
        
        % Create list of lipid species in the samples under standard conditions that were detected in D3:
        %includedInRatios(count,1) = lipid_i;
        includedControl(count) = i;
        
        %Verify inclusion of SD:
        getSD = absoluteLipids{i,4};
        if isempty(getSD)
            calculateSD = str2double(string(absoluteLipids(i,3)))/8;
            absoluteLipids{i,4} = cellstr(string(calculateSD));
        end
   
        % Obtain fold change and adjust absolute values for accessions
        for ii = 2:size(lipidRatios,1)
            for iii = 1:numel(col_i)
                colIdx = col_i(1,iii);
                % Calculate fold change for ecotype:
                FC = str2double(string(lipidRatios(ii,colIdx)));
                
                % Adjust abundance taking into account fold change
                ecotAbsoluteQ = controlAbsoluteQ * FC;

                % Verify that adjusted value don't exceed maximum lipid
                % content according to leaf distribution of dry weight 
                % reported by Li-Beisson et al. (as percent of total). Dry
                % weigth approx. 8% of fresh weigth (DOI:
                % 10.1199/tab.0133). In case that a value exceed the limit
                % of 0.008 g g-1 DW, the value that is assigned will be
                % equal to the conten under standard conditions.
                if ecotAbsoluteQ > 0.008
                    ecotAbsoluteQ = controlAbsoluteQ; 
                end
                if isempty(absoluteLipids{i,ii+3})
                    % Allocate abundance value if cell is empty 
                    absoluteLipids(i,ii+3) = cellstr(string(ecotAbsoluteQ));
                else
                    % If there are several isomers for the same species,
                    % sum the abundances to form a unique pool:
                    Qvalue_i = str2double(string(absoluteLipids{i,ii+3}));
                    absoluteLipids(i,ii+3) = cellstr(string(Qvalue_i + ecotAbsoluteQ));
                end
            end
        end
    end
end

% ###
% ##### 2). Add lipids with constant values to the list of absolute lipids:
% ###

% Note: There is a set of structural lipids which are part of the species 
% listed for samples under standard conditions. Those species that were not
% measured in the plants under darkness, are assumed to bear a constant
% value!

% Eliminate lipids measured under darkness from data set of constant lipids:
includedControl(includedControl == 0) = ''; 
constantLipids = QLipidsControl;
constantLipids(includedControl,:) = '';

% Eliminate lipid classes that are not part of the data set of constant lipids:
notConstant = {'PC', 'PE', 'PG', 'PS', 'PI', 'DG', 'TG', 'MGDG', 'DGDG', 'SQDG'};
constantAbb = extractBefore(constantLipids(:,1), ' ');
idxNonStrc = ismember(constantAbb, notConstant);
constantLipids(idxNonStrc == 1,:) = '';

% Add constant values to list of absolute lipids:
[idxAbsolute,idxConstant] = ismember(absoluteLipids(:,1), constantLipids(:,1));
idxAbsolute(1) = 0;
idxConstant(1) = 0;
constantS = constantLipids((idxConstant(idxConstant ~= 0)),3);
absoluteLipids((idxAbsolute==1),5:end) = repmat(constantS,1,((size(absoluteLipids,2)-4)));

% ###
% ##### 3). Find lipid species that were identified under standard conditions 
% ###       but were not detected in the samples exposed to stress:
missingIdx = cellfun(@isempty, absoluteLipids(:,5:end));
missingIdx = sum(missingIdx,2);
missingStress = absoluteLipids(missingIdx ~= 0, 1);
missingAbsolute = LipidsUnderDetection(missingStress, absoluteLipids, 1, lipidRatios);
absoluteLipids(missingIdx ~= 0, 5:end) = cellstr(string(missingAbsolute));

% ###
% ##### 4). Find lipid species that were identified in samples exposed to stress
% ###       but were not detected in the samples under standard conditions:
includedInRatios = QLipidsControl(includedControl,1);
missingControl = setdiff(lipidRatios(1,2:end),includedInRatios);
missingAbsolute = LipidsUnderDetection(missingControl', absoluteLipids, 2, lipidRatios);
missingControl = [missingControl',cell(numel(missingControl),1)];
missingControl = [missingControl, cellstr(string(missingAbsolute))];
absoluteLipids = [absoluteLipids;missingControl];

end
