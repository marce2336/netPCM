function referenceFluxes = getPoolSD(mt_i, referenceFluxes, processedData)
%##########################################################################
%
% This function compute the average standard deviation for a selected pool
% of lipids
%
%##########################################################################

% get list of lipid pools:
poolNames = fieldnames(referenceFluxes);

% get index for measurements specific to selected mutant:
idx_mti = strcmp(processedData.processedProfiles.SD_FC(1,:), mt_i);


for i = 1:numel(poolNames)
    pool_i = poolNames{i};
    
    switch pool_i
        case {'LPC','SPB','M_CoA','FA_CoA', 'FA_ACP','tCer','dCer','dCer18_1', 'dSP'}
            poolID = '] LPC ';
            
        case {'TG', 'FA'}
            poolID = '] TG ';
            
        case {'DG', 'MG', 'CDP_DG'}
            poolID = '] DG ';
            
        case {'PC','NMEthP','DMEthP','ChoP'}
            poolID = '] PC ';
            
        case {'LPA', 'PA', 'PGP','alkyl_ferulate'}
            % add here exception for mutant 'AT1G06520' for which data is
            % not available for LPA or PA:
            switch mt_i
                case {'AT1G06520'}
                    poolID = '] LPC ';
                otherwise
                    poolID = '] PA ';
            end
            
        case {'PE'}
            poolID = '] PE ';
            
        case {'SQNV'}
            poolID = '] SQDG ';
            
        case {'MGDG'}
            poolID = '] MGDG ';
            
        case {'DGDG'}
            poolID = '] DGDG ';
            
        case {'PG','G3P'}
            poolID = '] PG ';
    end
    
    % get standard deviation for measurements of selected species:
    idxSpecies = contains(processedData.processedProfiles.SD_FC(:,1), poolID);
    sd_mti = str2double(string(processedData.processedProfiles.SD_FC(idxSpecies,idx_mti)));

    % exclude SD values that are extremely high (>0.8)
    sd_mti = sd_mti(sd_mti<=0.8);
    
    % if all SD values are too high, take averages for all MTs:
    if isempty(sd_mti) 
        sd_mti = str2double(string(processedData.processedProfiles.SD_FC(idxSpecies,2:end)));
        sd_mti  = sd_mti (:);
        sd_mti = sd_mti(sd_mti<=0.8);
    end
    sd_mti = mean(sd_mti);
    
    % here we are assuming that the abundance of the lipid pool did not
    % change upon the mutation. Hence FC == 1.
    averageFC = 1;

    % add compyted sd and FC to reference fluxes data:
    referenceFluxes.(pool_i).sd_mti = sd_mti;
    referenceFluxes.(pool_i).averageFC = averageFC;
end



