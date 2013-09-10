function [conc_ss,v_ss,Vnet_ss] = solve_ss(expanded_model,y0,kinetic_param,perturbation)

% This function solves the ODEs and finds the rate of elementary and overall rxns
%
% INPUTS:
% -------
%             y0: Initial values of nomralized metabolite concentrations and enzyme fracitons
%  kinetic_param: kinetic parameters of rxns in expanded_model.rxn_f_b
%   perturbation: A vector whose size is equal to the number of enzymes. 
%                 It contains a value of
%                 1    if the enzyme is not perturbed
%                 0    if the enzyme is knocked out
%                 > 1  if the enzyme is overexpressed
% 
% OUTPUTS:
% --------
%       conc_ss: Steady state normalized metabolite concentrations and enzyme fractions
%         v_ss:  Rate of elementary rxns at steady state
%      Vnet_ss:  Rate of overall rxns at steady state
%
% Ali R. Zomorrodi- Costas Maranas Lab, April 2012
%

S_f_b = expanded_model.S_f_b;

% Solve system of equations

options = optimset('display','off','MaxIter',1000,'TolFun',1e-9,'TolX',1e-9);
[conc_ss,fval,exitflag]=fsolve(@mass_balance_ss,y0,options,expanded_model,kinetic_param,perturbation);

if exitflag < 0 
    return
elseif exitflag ~= 1
   fprintf('\nThe exit falg for solve = %i\n',exitflag);
end

%---------- Compute reaction rates ----------------
v_ss = zeros(size(kinetic_param));

% First compute the rate of reactions 
for j = 1:length(kinetic_param)
  v_ss(j) = kinetic_param(j);
  for i=1:size(S_f_b,1)
      if S_f_b(i,j) < 0
          v_ss(j)=v_ss(j)*conc_ss(i);
      end
  end
end

%---------- Compute Vnet ----------------
Vnet_ss = zeros(size(expanded_model.unexp_rxn_info,1),1);

for j=1:size(expanded_model.unexp_rxn_info,1)
   % If this is not an export reaction which is irreversible
   if expanded_model.unexp_rxn_info{j,6} ~= 3 
      % index of the first elementary rxn corresponding to j
      elem_index_first = expanded_model.unexp_rxn_info{j,2};

      % Corresponding indices of v_f and v_b. The first element is index of the forward
      % rxn whereas the second element is the index of the backward rxn
      v_indices=find(expanded_model.rxn_f_bInd_rxnInd==elem_index_first);

      % Vnet = v_f - v_b
      Vnet_ss(j) = v_ss(v_indices(1)) - v_ss(v_indices(2));

      % Check if v_f - v_b are equal for elemetnary steps
      for k=expanded_model.unexp_rxn_info{j,2}:expanded_model.unexp_rxn_info{j,3}
         v_k_indices=find(expanded_model.rxn_f_bInd_rxnInd==k);
         if abs(v_ss(v_k_indices(1)) - v_ss(v_k_indices(2)) - Vnet_ss(j)) > 0.1
%               fprintf('(v_f - v_b) not equal for all elementary steps of reaction %s at time %d\n',expanded_model.unexp_rxn_info{j,1},t);
         end
      end

   % if this an export rxn which is irreversilbe
   else
      % Vnet = v
      elem_index = expanded_model.unexp_rxn_info{j,2};
      v_index=find(expanded_model.rxn_f_bInd_rxnInd==elem_index);
      Vnet_ss(j) = v_ss(v_index);
   end
end

