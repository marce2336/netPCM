function ratios = AdjustMetContent(abundances, flag, commonMets)

switch flag
    case 1
        ratios = cell(size(abundances,1),size(abundances,2));
        ratios(:,1) = abundances(:,1);
        
        for i = 2:size(abundances,2)
            met = abundances(1,i);
            isMet = find(strcmp((string(commonMets(1,:))),(string(met))));
            
            if ~isempty(isMet)
                ratios(1,i) = met;
                
                for ii = 2:size(abundances,1)
                    accession = abundances(ii,1);
                    abundanceD0 = str2double(string(abundances(ii,i)));
                    isAccession = find(strcmp((string(commonMets(:,1))),(string(accession))));
                    
                    if ~isempty(isAccession)
                        abundanceD_i = str2double(string(commonMets(isAccession,isMet)));
                        ratio_i = cellstr(string(abundanceD_i/abundanceD0));
                        ratios(ii,i) = ratio_i;
                    end
                end
            end
        end

    case 2
        nRows           = size(abundances,1);
        ratios          = cell(nRows-1,size(abundances,2));
        ratios(2:end,1) = {'6h/Control', '21h/Control'};
        ratios(1,:)     = abundances(1,:);

        for i = 2:size(abundances,2)
            metControl_i = str2double(string(abundances(2,i)));
            met6h_i      = str2double(string(abundances(3,i)));
            met21h_i     = str2double(string(abundances(4,i)));
            ratio_6h     =  cellstr(string(met6h_i/metControl_i));
            ratio_21h    =  cellstr(string(met21h_i/metControl_i));
            ratios(2,i)  = ratio_6h;
            ratios(3,i)  = ratio_21h;
        end
end
end
