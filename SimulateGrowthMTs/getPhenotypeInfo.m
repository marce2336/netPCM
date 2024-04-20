function [phenoInfo, listKOs] = getPhenotypeInfo(flagMTset)
%##########################################################################
%
% This function use information from three databases (i.e. tair, SeedGenes,
% and ARALIPID), to search for phenotype description of the mutants. The
% results is a list with the phenotype description for each mutant.
%
% flagDataset  -  ('all') load information for the complete set of mutants
%                 ('filter') load information only for mutants corresponding 
%                            to genes included in the PLM
%                 ('lusk') load information for Arabidopsis mutants published by Lusk et al.
%                          (doi:https://doi.org/10.1093/pcp/pcac088)
%
%##########################################################################


% load list KOs:
pathFile = fullfile('InputDataSimulateGrowth','Phenotype_selected_KO_genes.xlsx');
switch flagMTset
    case {'filter'}
        listKOs = cellstr(string(readcell(pathFile, 'Sheet', 'Selected_KO_genes')));

    case {'all'}
        listKOs = cellstr(string(readcell(pathFile, 'Sheet', 'mutant_labels_and_descriptions')));
        
    case {'lusk'}
        listKOs = cellstr(string(readcell(pathFile, 'Sheet', 'lusk_mutants_set')));
end

% load tair info:
phenoTair = cellstr(string(readcell(pathFile, 'Sheet', 'Tair')));

% load SeedGenes info:
phenoSeedGenes = cellstr(string(readcell(pathFile, 'Sheet', 'SeedGenes')));
phenoSeedGenes(:,1) = strrep(phenoSeedGenes(:,1), 't', 'T');
phenoSeedGenes(:,1) = strrep(phenoSeedGenes(:,1), 'g', 'G');

% load ARALIPID info:
phenoAralipid = cellstr(string(readcell(pathFile, 'Sheet', 'ARALIPmutantDB_edited')));
phenoAralipid(:,1) = strrep(phenoAralipid(:,1), 't', 'T');
phenoAralipid(:,1) = strrep(phenoAralipid(:,1), 'g', 'G');


% 3). Identify phenotype info available for each KO:
phenoInfo = repmat({'NA'}, size(listKOs,1), 4);
phenoInfo(:,1) = listKOs(:,2);
phenoInfo(1,:) = {'locusID', 'tair', 'seedGenes', 'aralipid'};

phenoInfo = matchPhenotype(phenoInfo, listKOs(:,2), phenoTair, 2);
phenoInfo = matchPhenotype(phenoInfo, listKOs(:,2), phenoSeedGenes, 3);
phenoInfo = matchPhenotype(phenoInfo, listKOs(:,2), phenoAralipid, 4);


end