function newFormat = reformatModel(model,isPTS,cof_rxn_index,model_name);

% This function converts the model format provided by Jimmy the format used by codes written
% by Ali R. Zomorrodi in Costas' lab
%
% INPUT:
%-------
%          model: A matlab structure proivded by Dr. Liao's group containing 
%                 the following fields
%                          S: Stoichiomteric matrix 
%                       Sreg: A matrix representing regulation
%                  Vin_index: Index of the uptake reactions
%                 Vout_index: Index of the export reactions
%                       SGFE: Standard free Gibbs energy of each reaction
%                    EnzName: Enzyme names
%                  MetabName: Metabolite names
%                   GFERange: Range of the free Gibbs energy
%                 MetabRange: Range of the metabolite concentrations
%                 Vcof_index: Index oc the cofators
%                        Sbm: Stoichiomteric coefficient of the biomass
%                   rVnet_bm: Net rate of the biomass (e.g., for ATP and ADP we need 
%                             the index of ATP)
%                      rVnet: Net rate of fluxes at steady state
%              Perturbations: A matrix representing perturbations in the model. 
%                             Rows correspond to to genes and columns to each perturbation
%
%
%  cof_rxn_index: A vector containing the index of rxns responsible for 
%                 cofactor interconversions (e.g., 
%                 ATP <--> ADP, NADH <--> NAD,   NADPH <--> NADP)
%          isPTS: 0 if there is no pts rxn in the model and a vector with two elements
%                 where the first element takes a value of one representin there is a 
%                 pts rxn in the model and the second element is the index of pts rxn in 
%                 input data structure
%     model_name: A string containing the name of the model
%
% OUTPUTS:
% ---------
%          metab: Name of the metabolites
%            rxn: Name of the rxns 
%              S: Stoichiomteric matrix        
%         Enzyme: Name of enzymes
%      enzymeReg: A matrix represenitng which metabolite affect the activity of which
%                 enzyme. Rows correspond to metabolites and columns correspon to enzymes
%                 If element (i,j) of this matrix is non-zeor it implies that metabolite
%                 if has an inhibiton/activation effect on enzyme j:
%                 -1 = competetive inhibution  -2 = Uncompetetive inhibition
%                 -3 = mixed inhibution        -4 = activation 
%        rxntype: Type of rxn: 1 = internal  2 = uptake 3 = export 4 = PTS 
%                 5 = cofactor convert (e.g., ATP <--> ADP)
%           Vss:  Steady state fluxes. The first column correspond to the wild-type.
%                 The rest of the columns correspond to perturbations 
%           GPR:  A matrix representing the gene-(protein)-reaction associatons.
%                 Rows correspond to genes and columns to reactions
%          SGFE:  Standard Gibbs free energy
%    metabRange:  Range of metabolite concentrations                  
%     cof_index:  Index oc the cofators
% perturb_names:  Name of the perturbed strain. The first element is wild showing the
%                 wild type strain (i.e., no perturbations)     
% Perturbations:  A matrix representing perturbations in the model. Rows correspond to
%                to genes and columns to each perturbation
%

% Metabolite names
newFormat.metab = model.MetabName;

% Reaction names
newFormat.rxn = model.EnzName;

% Stoichiomteric matrx
newFormat.S=model.S;

% Name of enzymes
newFormat.Enzyme = model.EnzName;

% Gene (enzyme) reaction association
newFormat.GPR=zeros(max(size(newFormat.Enzyme)),max(size(newFormat.rxn)));

for i=1:size(newFormat.GPR,1)
  for j=1:size(newFormat.GPR,2)
     if (i==j) && isempty(find(model.Vout_index == j))
%     if (i==j) && (j ~= model.Vout_index)
           newFormat.GPR(i,j)=1;
     end
  end
end

% Enzyme regulation
% -1 = competitive inhibition   -2 = uncompetitive inhibition
% -3 = mixed inhibition         -4 = activation
newFormat.enzymeReg =zeros(max(size(newFormat.metab)),max(size(newFormat.Enzyme))); 

for j = 1:size(newFormat.enzymeReg,2)
  for i = 1:size(newFormat.enzymeReg,1)
        if model.Sreg(i,j) ~= 0
           newFormat.enzymeReg(i,j) = model.Sreg(i,j);
        end
  end
end


% rxntype
newFormat.rxntype=ones(size(newFormat.rxn));
newFormat.rxntype(model.Vin_index)=2;
newFormat.rxntype(model.Vout_index)=3;


% Steady state flux
newFormat.Vss = model.rVnet;

% Standard free Gibbs energy 
newFormat.SGFE = model.SGFE;

% Range of metabolite concentrations
newFormat.metabRange = model.MetabRange;

% Coefficient of the primary cofactors (in the field for metabolite names)
newFormat.cof_metab_index = model.Vcof_index;

% Coefficineit of rxns converting cofactors (e.g., ATP <--> ADP)
newFormat.cof_rxn_index = cof_rxn_index; 

%--------- Perturbations -------------
newFormat.Perturbations = ones(max(size(newFormat.Enzyme)),1); 


% Find which column of model.Perturbations correspond to the wildtype and remove them
counter = 0;

wild_index=[];%+++++++++++++
for j=1:size(model.Perturbations,2)
   % If there is no perturbation in column j
   if isempty(find(model.Perturbations(:,j) ~= 1))
       counter = counter + 1;
       if counter == 1
         wild_index = j;
       %end%++++++++++++
       else
         wild_index = [wild_index j] 
       end
   end
end

model.Perturbations(:,wild_index) = [];

if max(size(model.Perturbations)) > 0
  
newFormat.Perturbations(model.Vin_index(1):model.Vout_index(end),2:size(model.Perturbations,2)+1) = 1; 
   newFormat.Perturbations(1:model.Vin_index(1) - 1,2:size(model.Perturbations,2)+1) = model.Perturbations; 
end

% Perturbation names
newFormat.perturb_names = createPerturbName(newFormat); 

% Implement PTS system
if isPTS(1)==1
  enzyme_num = length(model.EnzName);
  pts_rxn_index = isPTS(2);
  newFormat.rxntype(pts_rxn_index)=4;

  % Replace the name pts with EI
  newFormat.Enzyme(pts_rxn_index)={'EI'};
 
  newFormat.Enzyme(enzyme_num+1)={'HPr'};
  newFormat.GPR(enzyme_num+1,pts_rxn_index)=1;
  newFormat.Perturbations(enzyme_num+1,:)=1;

  newFormat.Enzyme(enzyme_num+2)={'EIIA'};
  newFormat.GPR(enzyme_num+2,pts_rxn_index)=1;
  newFormat.Perturbations(enzyme_num+2,:)=1;

  newFormat.Enzyme(enzyme_num+3)={'EIIBC'};
  newFormat.GPR(enzyme_num+3,pts_rxn_index)=1;
  newFormat.Perturbations(enzyme_num+3,:)=1;
  
end

newFormat.model_name=model_name;