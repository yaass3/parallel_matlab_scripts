function [T,conc_perturb,v_perturb,Vnet_perturb] = perturbModel(expanded_model,kinetic_param,t_interval,sampled_e_vec,np)

% This mode perturbs a model in the ensemble
%
% INPUTS:
% -------
%       expanded_model: Expanded model
%                   np: The column number of the filed Perturbation in expanded model
%
% OUTPUTS:
% --------
%
%
% Ali R. Zomorrodi, Costas Maranas Lab, April 2012
%

% Set the initial concentrations
y0=sampled_e_vec;
y0(expanded_model.cof_metab_index) = 1;
y0(expanded_model.cof_metab_index+1) = 1;

% Other metabolites in the unexpanded model that are not cofactors
for i = 1:cell2mat(expanded_model.unexp_rxn_info(1,2))
  if isempty(find(expanded_model.cof_metab_index == i)) && isempty(find(expanded_model.cof_metab_index+1 == i))
     % Use random numbers to do the plots
     % y0(1:expanded_model.cof_metab_index(1) -1) = rand;   

     % Use 1 otherwise 
     y0(1:expanded_model.cof_metab_index(1) -1) = 1;
  end
end

% Index of the perturbed enzymes
perturb_enz_indices = find(expanded_model.Perturbations(:,np) ~= 1); 
for ni=1:length(perturb_enz_indices)
    % Find the index of the free enzyme
    free_enz_index = expanded_model.metabIndexInfo{1,3} + perturb_enz_indices(ni);

    % Find the indices of corresponding enzyme complex
    enzComplex_indices = find(expanded_model.enz_enzComplex(:,perturb_enz_indices(ni)) ~= 0);

    %--- knockout ---
    if expanded_model.Perturbations(perturb_enz_indices(ni),np) == 0

       % Initial condition for this enzyme and its corresponding enzyme complexes are zero
       y0(free_enz_index) = 0;
       y0(enzComplex_indices) = 0;              

    %--- overexpression ---
    elseif expanded_model.Perturbations(perturb_enz_indices(ni),np) > 1

       % Sum of the enzyme fractions should be equal to the value of perturbation which
       % is greater than. This sum is currently one. So, if add the difference of the value
       % of perturbaiton and one to the initial value of the enzyme fraction for free enzyme
       % the sum will be equal to the value of perturbation 
       y0(free_enz_index) = y0(free_enz_index) + [expanded_model.Perturbations(perturb_enz_indices(ni),np) - 1];

    % if downregulaiton
    elseif expanded_model.Perturbations(perturb_enz_indices(ni),np) < 1
       % Sample e between zero and the perturbation value 
       % Total number of enzyme fractions associated with each enzyme
       enz_frac_num = 1 + length(enzComplex_indices); 

       % Perform the sampling
       e_sample_pool=rand(10000,enz_frac_num);

       % Normalize to the value of perturbation
       e_sample_pool = expanded_model.Perturbations(perturb_enz_indices(ni),np)*e_sample_pool./repmat(sum(e_sample_pool,2),1,enz_frac_num);

       % Draw one sample from the sample pool randomly
       isEqualToPerturbValue = 0;
       while isEqualToPerturbValue == 0
       random_index=randi(10000,1);
         if sum(e_sample_pool(random_index,:)) == expanded_model.Perturbations(perturb_enz_indices(ni),np)
           isEqualToPerturbValue = 1;
         end
       end  % end of while
       y0(free_enz_index) = e_sample_pool(random_index,1);
       y0(enzComplex_indices) = e_sample_pool(random_index,2:end);

    end % end of if expanded_model.Perturbations(perturb_enz_indices(ni),np) == 0 else ...
end   % end of for ni

%---- Now solve the ODEs ----
% Solve the system of ODEs
[T,conc_perturb,v_perturb,Vnet_perturb] = solve_ode(expanded_model,t_interval,y0,kinetic_param);

