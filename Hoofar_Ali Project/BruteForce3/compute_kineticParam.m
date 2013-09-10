function  kinetic_param = compute_kineticParam(expanded_model,sampled_e_vec,elem_rate)

% This function computes the kinetic parameters 
%
% INPUTS:
% -------
% expanded_model: Expanded model
%  sampled_e_vec: A vector whose size is equal to the number of metabolites in 
%                 expanded network and contains the enzyme fracitons 
%      elem_rate: Rate of elementary reactions
%
% OUTPUTS:
% --------
% kinetic_param: Kinetic parameters for each reaciton in rxn_f_b
%


kinetic_param = zeros(size(expanded_model.rxn_f_b));

for i=1:length(expanded_model.rxn_f_b)
  % if this is an elementary (and hence enzyme-catalyzed) rxn
  %  v = k*e  or v = k*e1*e2*...
  if expanded_model.isElementary(expanded_model.rxn_f_bInd_rxnInd(i)) == 1
     % Indices of reactants of the this reaction
     reactant_indices = find(expanded_model.S_f_b(:,i) < 0 );
     enz_enzComplex_indices = reactant_indices(find(reactant_indices >= expanded_model.metabIndexInfo{2,2}));
     kinetic_param(i) = elem_rate(i)/prod(sampled_e_vec(enz_enzComplex_indices));
  % If not an elementary (and hence non-enzymatic) rxn
  % v = k
  else
     kinetic_param(i) = elem_rate(i);
  end 
end

