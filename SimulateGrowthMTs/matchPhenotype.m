function phenoInfo = matchPhenotype(phenoInfo, KOs, database, col)
%##########################################################################
%
% This function compares the list of Arabidopsis KOs againts a list of
% mutants with phenotype description from different database resources.
% When a match is found, the corresponding index of the description is
% retrieved.
%
%##########################################################################


[isPhenoGene,locb] = ismember(KOs, database(:,1));

locb(locb == 0) = '';

switch col
    case 2
        phenoInfo(isPhenoGene,col) = database(locb,2);

    case 3
        phenoInfo(isPhenoGene,col) = database(locb,10);

    case 4
        phenoInfo(isPhenoGene,col) = database(locb,12);
end

end