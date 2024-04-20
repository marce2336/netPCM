function CalculateRatios

%% Calculation of ratios of metabolites for dark- vs. light-grown Arabidopsis accessions:
%  This script allows the calculation of the ratios among dark- vs. 
%  light-grown Arabidopsis accessions for relative abundance data obtained
%  for soluble intracellular metabolites and lipids.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% 1). Load relative abundance data sets:
    pathMets = fullfile('InputData','NormalizedIntensities.xlsx');
    [~, ~, Day0] = xlsread(pathMets,'0d','A1:HP285');
    Day0 = RemoveIsNaN(Day0);

    [~, ~, Day3] = xlsread(pathMets,'3d','A1:HF285');
    Day3 = RemoveIsNaN(Day3);

    [~, ~, Day6] = xlsread(pathMets,'6d','A1:HE288');
    Day6 = RemoveIsNaN(Day6);

    [~, ~, Col0] = xlsread(pathMets,'Col-0','A1:FF4');
    Col0 = RemoveIsNaN(Col0);

    %% 2). Find common metabolites detected among D0 and dark-treated samples:
    CommonMetsDay3 = cell(size(Day3,1),size(Day3,2)); %List of common mets among Day0 and Day3
    CommonMetsDay3(:,1) = Day3(:,1);

    CommonMetsDay6 = cell(size(Day6,1),size(Day6,2)); %List of common mets among Day0 and Day6
    CommonMetsDay6(:,1) = Day6(:,1);

    DataSetsList = {Day3,Day6};
    Count1 = 1;
    Count2 = 1;

    for i = 1:length(DataSetsList)
        DataSet = DataSetsList{1,i};

        for ii = 2:size(DataSet,2)
            Met = DataSet(1,ii);
            IsMet = find(strcmp((string(Day0(1,:))),(string(Met))));
            if ~isempty(IsMet) && i == 1
                Count1 = Count1 + 1;
                CommonMetsDay3(1,Count1) = Met;
                CommonMetsDay3(2:(size(Day3,1)),Count1) =  Day3(2:size(Day3,1),ii);

            elseif ~isempty(IsMet) && i == 2
                Count2 = Count2 + 1;
                CommonMetsDay6(1,Count2) = Met;
                CommonMetsDay6(2:(size(Day6,1)),Count2) =  Day6(2:size(Day6,1),ii);
            end

        end   
    end

    clearvars i ii Met IsMet Count1 Count2 DataSetsList DataSet ans

    CommonMetsDay3(:,all(cellfun(@isempty, CommonMetsDay3),1)) = []; 
    CommonMetsDay6(:,all(cellfun(@isempty, CommonMetsDay6),1)) = [];


    %% 3). Find metabolites that differ among Day0 and Day3 samples:

    MissingD3 = setdiff(Day0(1,:),Day3(1,:)); %List of metabolites present in D0 but not in D3
    MissingD0 = setdiff(Day3(1,:),Day0(1,:)); %List of metabolites present in D3 but not in D0

    % ADDING MISSING METABOLITES TO D0  REGARDING D3 SAMPLES:
    % The metabolites that are present in D3 but not in D0 samples, correspond
    % to metabolites accumulated under the corresponding stress condition but
    % were absent or at very low levels under the control condition, since they
    % were nonexistent or were probably under the detection limit. Therefore, 
    % the missing metabolites will be added to the D0 list with an abundance
    % magnitude equal to 1/5 of the minimum positive abundance value detected 
    % in all the metabolites among all the accessions evaluated, corresponding 
    % to 1.9378e-05.  

    Day0_RatiosD3 = [Day0,cell((size(Day0,1)),(length(MissingD0)))];
    Day0_RatiosD3(1,(size(Day0,2)+1):size(Day0_RatiosD3,2)) = MissingD0;
    Day0_RatiosD3(2:end,(size(Day0,2)+1):size(Day0_RatiosD3,2)) = {'1.9378e-05'};

    % The remaining missing values are replaced by 1/5 of the minimum positive
    % abundance value of their corresponding variables.

    IsEmpty = cellfun(@isempty, Day0_RatiosD3);

    for i = 2:size(Day0_RatiosD3,2)
        Variable = str2double(string(Day0_RatiosD3(2:end,i)));
        MinValue = (min(Variable))/5;
        MissingValues = IsEmpty(:,i) == 1;
        Day0_RatiosD3((MissingValues == 1),i) = cellstr(string(MinValue));
    end

    %The same metabolites will be added to the CommonMetsDay3 list with the
    % corresponding abundance values.
    OldLength = size(CommonMetsDay3,2);
    CommonMetsDay3 = [CommonMetsDay3,cell(size(CommonMetsDay3,1),length(MissingD0))];

    for i = 1:length(MissingD0)
        Met = MissingD0(1,i);
        IsMet = find(strcmp((string(Day3(1,:))),(string(Met))));
        if ~isempty(IsMet)
            CommonMetsDay3(:,OldLength+i) = Day3(:,IsMet);
        end
    end

    % ADDING MISSING METABOLITES TO D3 REGARDING D0 SAMPLES:
    % The metabolites that are present in D0 but not in D3 samples, correspond
    % to metabolites accumulated in the control condition but probably degraded
    % under the corresponding stress condition, since they were nonexistent 
    % or were probably under the detection limit. Therefore,the missing
    % metabolites will be added to the D3 list with an abundance magnitude equal
    % to 1/5 of the minimum positive abundance value detected in all the
    % metabolites among all the accessions evaluated, corresponding to 1.9794e-06.

    CommonMetsDay3 = [CommonMetsDay3,cell((size(CommonMetsDay3,1)),(length(MissingD3)))];
    CommonMetsDay3(1,(size(CommonMetsDay3,2)-length(MissingD3)+1):end) = MissingD3;
    CommonMetsDay3(2:end,(size(CommonMetsDay3,2)-length(MissingD3)+1):end) = {'1.9794e-06'};

    % The remaining missing values are replaced by 1/5 of the minimum positive
    % abundance value of their corresponding variables.

    IsEmpty = cellfun(@isempty, CommonMetsDay3);

    for i = 2:size(CommonMetsDay3,2)
        Variable = str2double(string(CommonMetsDay3(2:end,i)));
        MinValue = (min(Variable))/5;
        MissingValues = IsEmpty(:,i) == 1;
        CommonMetsDay3((MissingValues == 1),i) = cellstr(string(MinValue));
    end

    %% 4). Find metabolites that differ among Day0 and Day6 samples:

    MissingD6 = setdiff(Day0(1,:),Day6(1,:)); %List of metabolites present in D0 but not in D6
    MissingD0 = setdiff(Day6(1,:),Day0(1,:)); %List of metabolites present in D6 but not in D0

    % The metabolites that are present in D6 but not in D0 samples, correspond
    % to metabolites accumulated under the corresponding stress condition but
    % were absent or at very low levels under the control condition, since they
    % were nonexistent or were probably under the detection limit. Therefore, 
    % the missing metabolites will be added to the D0 list with an abundance
    % magnitude equal to 1/5 of the minimum positive abundance value detected 
    % in all the metabolites among all the accessions evaluated, corresponding
    % to 1.9378e-05.

    Day0_RatiosD6 = [Day0,cell((size(Day0,1)),(length(MissingD0)))];
    Day0_RatiosD6(1,(size(Day0,2)+1):size(Day0_RatiosD6,2)) = MissingD0;
    Day0_RatiosD6(2:end,(size(Day0,2)+1):size(Day0_RatiosD6,2)) = {'1.9378e-05'};

    % The remaining missing values are replaced by 1/5 of the minimum positive
    % abundance value of their corresponding variables.

    IsEmpty = cellfun(@isempty, Day0_RatiosD6);

    for i = 2:size(Day0_RatiosD6,2)
        Variable = str2double(string(Day0_RatiosD6(2:end,i)));
        MinValue = (min(Variable))/5;
        MissingValues = IsEmpty(:,i) == 1;
        Day0_RatiosD6((MissingValues == 1),i) = cellstr(string(MinValue));
    end

    %The same metabolites will be added to the CommonMetsDay6 list with the
    % corresponding abundance values.
    OldLength = size(CommonMetsDay6,2);
    CommonMetsDay6 = [CommonMetsDay6,cell(size(CommonMetsDay6,1),length(MissingD0))];

    for i = 1:length(MissingD0)
        Met = MissingD0(1,i);
        IsMet = find(strcmp((string(Day6(1,:))),(string(Met))));
        if ~isempty(IsMet)
            CommonMetsDay6(:,OldLength+i) = Day6(:,IsMet);
        end
    end

    % ADDING MISSING METABOLITES TO D6 REGARDING D0 SAMPLES:
    % The metabolites that are present in D0 but not in D6 samples, correspond
    % to metabolites accumulated in the control condition but probably degraded
    % under the corresponding stress condition, since they were nonexistent 
    % or were probably under the detection limit. Therefore,the missing
    % metabolites will be added to the D6 list with an abundance magnitude equal
    % to 1/5 of the minimum positive abundance value detected in all the
    % metabolites among all the accessions evaluated, corresponding to 1.193e-10.

    CommonMetsDay6 = [CommonMetsDay6,cell((size(CommonMetsDay6,1)),(length(MissingD6)))];
    CommonMetsDay6(1,(size(CommonMetsDay6,2)-length(MissingD6)+1):end) = MissingD6;
    CommonMetsDay6(2:end,(size(CommonMetsDay6,2)-length(MissingD6)+1):end) = {'1.193e-10'};

    % The remaining missing values are replaced by 1/5 of the minimum positive
    % abundance value of their corresponding variables.

    IsEmpty = cellfun(@isempty, CommonMetsDay6);

    for i = 2:size(CommonMetsDay6,2)
        Variable = str2double(string(CommonMetsDay6(2:end,i)));
        MinValue = (min(Variable))/5;
        MissingValues = IsEmpty(:,i) == 1;
        CommonMetsDay6((MissingValues == 1),i) = cellstr(string(MinValue));
    end

    clearvars Variable MinValue MissingValues OldLength IsEmpty

    %% 5). Calculate ratios:
    % 
    % Calculate Day3 vs. Control(day0) ratios:
    RatiosD3 = AdjustMetContent(Day0_RatiosD3, 1, CommonMetsDay3);

    % Calculate Day6 vs. Control(day0) ratios:
    RatiosD6 = AdjustMetContent(Day0_RatiosD6, 1, CommonMetsDay6);

    % Calculate ratios for Col-0:
    RatiosCol0 = AdjustMetContent(Col0, 2);
    
    % Save output variables:
    pathMets = fullfile('OutputData','Ratios.xlsx');
    xlswrite(pathMets,RatiosD3,'RatiosD3')
    xlswrite(pathMets,RatiosD6,'RatiosD6')
    xlswrite(pathMets,RatiosCol0,'RatiosCol0')
    
end

