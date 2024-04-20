function [ChainCoeff, BackboneCoeff] = getAccessionCoefficients(ChainMets, BackboneMets, SLIMEcoeff)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'getAccessionCoefficients' function calculates the stoichiometric 
%  coefficients for the chains and backbones of the lipid pseudoreaction:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get chain coefficients of metabolites of lipid pseudoreactions in model:
ChainCoeff = zeros(size(ChainMets,1),1);

for i = 1:size(ChainMets,1)
     ChainID = ChainMets{i,1};
     isChain = strcmp(SLIMEcoeff.chainData.metNames,ChainID);
     
     if sum(isChain) == 1
         ChainCoeff(i,1) = SLIMEcoeff.chainData.abundance(isChain);
     else
         ChainCoeff(i,1) = 0;
     end
     
 end

% Get backbone coefficients of metabolites of lipid pseudoreactions in model:
BackboneCoeff = zeros(size(BackboneMets,1),1);
 for i = 1:size(BackboneMets,1)
     BackboneID = BackboneMets{i,1};
     isBackbone = strcmp(SLIMEcoeff.lipidData.metNames,BackboneID);
     
     if sum(isBackbone) == 1
         BackboneCoeff(i,1) = SLIMEcoeff.lipidData.abundance(isBackbone);
     else
         BackboneCoeff(i,1) = 0;
     end 
 end
end