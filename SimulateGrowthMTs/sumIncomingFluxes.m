function referenceFluxes = sumIncomingFluxes(referenceFluxes, flagFieldType)
%##########################################################################
%
% The flux sum of reactions (FSR) is computed by adding the incoming fluxes
% to each lipid pool (v_i^Col0). The incoming fluxes were in turn obtained 
% preliminary for Arabidopsis wild-type plants via pFBA:
% this the following formula is implemented:
% v_i^Col0  =  v_(rxn_1) + v_(rxn_2) + ... + v_(rxn_n)
%
% flagFieldType - (0). Default. The incoming fluxes are stored in nested struct fields
%                 (1) The incoming fluxes are not stored as nested fields
%
%##########################################################################

if nargin < 2
    flagFieldType = 0;
end

listLipids = fieldnames(referenceFluxes);


switch flagFieldType
    case 0
        for i = 1:numel(listLipids)
            lipid_i = listLipids{i};
            speciesPool_i = fieldnames(referenceFluxes.(lipid_i));
            idxLabel = strcmp(speciesPool_i,'speciesID');
            speciesPool_i(idxLabel) = '';
            fsr_i = 0;

            for j = 1:numel(speciesPool_i)
                species_ij = speciesPool_i{j};
                fluxesSpecies_ij = referenceFluxes.(lipid_i).(species_ij).incomingFluxes;
                fsr_i = fsr_i + sum(fluxesSpecies_ij);
            end

            referenceFluxes.(lipid_i).fsr = fsr_i;
        end

    case 1
        for i = 1:numel(listLipids)
            lipid_i = listLipids{i};
            fluxesSpecies_ij = referenceFluxes.(lipid_i).incomingFluxes;
            fsr_i = sum(fluxesSpecies_ij);
            referenceFluxes.(lipid_i).fsr = fsr_i;
        end
end
end


