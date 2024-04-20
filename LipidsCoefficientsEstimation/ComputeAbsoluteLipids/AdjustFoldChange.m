function missingAbsoluteNew = AdjustFoldChange(missingMets, missingAbsolute, lipidRatios, absoluteLipids)

for i = 1:numel(missingMets)
            missing_i = missingMets{i};
            controlAbsoluteQ = missingAbsolute(i,1);
            idxLipid_i = find(strcmp(lipidRatios(1,:), missing_i));

            for j = 2:size(lipidRatios,1)

                for ij = 1:numel(idxLipid_i)
                    isomer = idxLipid_i(ij);
                    ecotype_i = lipidRatios{j,1};
                    FC_i = str2double(string(lipidRatios(j,isomer)));
                    idxEcotype = (find(strcmp(absoluteLipids(1,:), ecotype_i)));
                    
                    ecotAbsoluteQ = controlAbsoluteQ * FC_i; % Adjust abundance taking into account fold change
                    
                    if ecotAbsoluteQ > 0.008
                        ecotAbsoluteQ = controlAbsoluteQ; 
                    end
    
                    if missingAbsolute(i,idxEcotype-2) == 0
                        % Allocate abundance value if cell is empty 
                        missingAbsolute(i,idxEcotype-2) = ecotAbsoluteQ; 
                    else
                        % If there are several isomers for the same species, sum up the abundances to form a unique pool
                        Qvalue_i = missingAbsolute(i,idxEcotype-2);
                        missingAbsolute(i,idxEcotype-2) = Qvalue_i + ecotAbsoluteQ;
                    end
                end
            end
end

missingAbsoluteNew = missingAbsolute;

end
