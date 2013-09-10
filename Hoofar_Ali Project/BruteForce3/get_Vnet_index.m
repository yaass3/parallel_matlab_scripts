function   unexp_rxn_index = get_Vnet_index(expanded_model,i)

% Function to find the index of corresponding rxn in unexpanded model 
%
%   expanded_mode: expanded model
%               i: Index of the elementary rxn in expanded_model.rxn
% unexp_rxn_index: Index of the corresponding rxn in unexpanded network
%

% start_indices_vec (which start indices are less than or equal to i)
start_ind_vec = cell2mat(expanded_model.unexp_rxn_info(:,2)) <= i;

% End (which end indices are greater than or equal to i)
end_ind_vec = cell2mat(expanded_model.unexp_rxn_info(:,3)) >= i;

% for which rxn in the unexpanded model both elements of start_ind_vec and
% start_ind_vec are equal (they are equal only if they are both one)
unexp_rxn_index = find((start_ind_vec == end_ind_vec) == 1);

