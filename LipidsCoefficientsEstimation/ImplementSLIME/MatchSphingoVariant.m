function [idxs, chains] = MatchSphingoVariant(sphingoVariant, matchedMets, LCB, shortTail, idxs, backCode)

splitshortTail = strsplit(shortTail{1,1},':');
LCBtail = strsplit(LCB{1,1},':');

% List with sphingolipids variants:
         % LCB-ID | #OH-groups | Full-structure | LCB-DB | #C-atoms-LCB | OH-shortTail | shortTail-ID   
sphingos = {'d'        '2'          '1OH,3OH'       '0'        '18'           '0'         ''        %sphingolipid-1: Cer d18:0-Cxx:x
            'd'        '2'          '1OH,3OH'       '1'        '18'           '0'         ''        %sphingolipid-2: Cer d18:1-Cxx:x
            'd'        '2'          '1OH,3OH'       '0'        '18'           '1'         'h'       %sphingolipid-3: Cer d18:0-hCxx:x
            'd'        '2'          '1OH,3OH'       '1'        '18'           '1'         'h'       %sphingolipid-4: Cer d18:1-hCxx:x
            'd'        '2'          '1OH,3OH'       '2'        '18'           '0'         ''        %sphingolipid-2: Cer d18:2-Cxx:x
            'd'        '2'          '1OH,3OH'       '2'        '18'           '1'         'h'       %sphingolipid-2: Cer d18:2-hCxx:x
            't'        '3'          '1OH,3OH,4OH'   '0'        '18'           '0'         ''        %sphingolipid-5: Cer t18:0-Cxx:x
            't'        '3'          '1OH,3OH,4OH'   '1'        '18'           '0'         ''        %sphingolipid-6: Cer t18:1(E)-Cxx:x
            't'        '3'          '1OH,3OH,4OH'   '0'        '18'           '1'         'h'       %sphingolipid-7: Cer t18:0-hCxx:x
            't'        '3'          '1OH,3OH,4OH'   '1'        '18'           '1'         'h'  };   %sphingolipid-8: Cer t18:1(E)-hCxx:x

switch sphingoVariant
    case 1
        % Find matching sphingolipid in model metabolites:
        IsLipid = strcmp(backCode,matchedMets(:,2)).*contains(matchedMets(:,3),LCB{1,1}).*contains(matchedMets(:,3),LCB{1,2})...
            .*contains(matchedMets(:,4),shortTail{1,1}) == 1;
        matchedMets = matchedMets(IsLipid, :);
        idxs = idxs(IsLipid, :);
        
        % Double check features of short tail:
        [idxs] = VerifyshortTail(shortTail, matchedMets, idxs);
        
        % Find features of chains for sphingolipid and match to the corresponding variant: 
        LCBid = strcmp(LCB{1,2},sphingos(:,3)).*strcmp(LCBtail{1,2},sphingos(:,4)).*strcmp('0',sphingos(:,6)) == 1;
        LCBchain = [sphingos{LCBid,1},LCBtail{1,1},':',LCBtail{1,2}];
        shortChain = [sphingos{LCBid,7},splitshortTail{1,1},':',splitshortTail{1,2}];
        chains = [cellstr(LCBchain),cellstr(shortChain)];
        
    case 2
        % Find matching sphingolipid in model metabolites:
        IsLipid = strcmp(backCode,matchedMets(:,2)).*contains(matchedMets(:,3),LCB{1,1}).*contains(matchedMets(:,3),LCB{1,2})...
            .*contains(matchedMets(:,4),shortTail{1,1}).*contains(matchedMets(:,4),shortTail{1,2}) == 1;
        matchedMets = matchedMets(IsLipid, :);
        idxs = idxs(IsLipid, :);
        
        % Double check features of short tail:
        [idxs] = VerifyshortTail(shortTail, matchedMets, idxs);
        
        % Find features of chains for sphingolipid and match to the corresponding variant:
        LCBid = strcmp(LCB{1,2},sphingos(:,3)).*strcmp(LCBtail{1,2},sphingos(:,4)).*strcmp('1',sphingos(:,6)) == 1;
        LCBchain = [sphingos{LCBid,1},LCBtail{1,1},':',LCBtail{1,2}];
        shortChain = [sphingos{LCBid,7},'C',splitshortTail{1,1},':',splitshortTail{1,2}];
        chains = [cellstr(LCBchain),cellstr(shortChain)];
        
    case 3
        % Find matching sphingolipid in model metabolites:
        IsLipid = strcmp(backCode,matchedMets(:,2)).*contains(matchedMets(:,3),LCB{1,1}).*contains(matchedMets(:,4),shortTail{1,1})...
            .*contains(matchedMets(:,4),shortTail{1,2}) == 1;
        matchedMets = matchedMets(IsLipid, :);
        idxs = idxs(IsLipid, :);
        
        % Double check features of short tail:
        [idxs] = VerifyshortTail(shortTail, matchedMets, idxs);
        
        % Find features of chains for sphingolipid and match to the corresponding variant:
        LCBid = strcmp(extractBefore(LCB{1,1},'18'),sphingos(:,1)).*strcmp(LCBtail{1,2},sphingos(:,4)).*strcmp('1',sphingos(:,6)) == 1;
        LCBchain = [LCBtail{1,1},':',LCBtail{1,2}];
        shortChain = [sphingos{LCBid,7},'C',splitshortTail{1,1},':',splitshortTail{1,2}];
        chains = [cellstr(LCBchain),cellstr(shortChain)];
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function [idxs] = VerifyshortTail(shortTail, matchedMets, idxs)

% Check if short tail contains functional groups and act accordingly:
functionalGroups = size(shortTail, 2);

switch functionalGroups
    % If short tail doesn't contain functional groups eliminate matching
    % metabolites with functional groups from list:
    case 1
        isFunctionalGroup = contains(matchedMets(:,4), 'OH');
        idxs = idxs(isFunctionalGroup == 0, :);
        
    % If short tail contains functional groups eliminate matching
    % metabolites without functional groups from list:
    case 2
        isFunctionalGroup = contains(matchedMets(:,4), 'OH');
        idxs = idxs(isFunctionalGroup == 1, :);
end

end