function [T,conc,v,Vnet] = solve_ode(expanded_model,t_interval,y0,kinetic_param)

% This function solves the ODEs and finds the rate of elementary and overall rxns
% t_interval: The time interval for simulations
%            y0: initial values of x and e
%                For cofactors y0 should be equal to one
%                For enzyme fractions use the sampled e values
% kinetic_param: kinetic parameters of rxns
% 
% Outputs:
%             T: Time vector
%          conc: Normalized metabolite concentrations and enzyme fractions:
%                Columns correspond to the time points 
%                Rows correspond to metabolites 
%             v: Rate of elementary rxns at each time point
%                Columns correspond to time points
%                Rows correspond to reactions
%          Vnet: Rate of overall rxns in the unexpanded model at each time point
%                Columns correspond to time points
%                Row correspond to reactions
%
% Ali R. Zomorrodi- April 2012
%

S_f_b = expanded_model.S_f_b;


% Solve system of ODEs
eval(strcat('options=odeset(''Jacobian'',@(t,y)jacobFun_',expanded_model.model_name,'(t,y,kinetic_param));'))
[T,Y]=ode15s(@(t,y)mass_balance_ode(t,y,kinetic_param,S_f_b),t_interval,y0,options);

%---------  Concentrations ---------------
% Rows in conc correspond to metabolites and columns correspond to time points
conc = Y'; 

%---------- Compute reaction rates ----------------
% Total number of time points
time_point_num=length(T);

v = zeros(length(kinetic_param),time_point_num);

% First compute the rate of reactions 
for j = 1:length(kinetic_param)
 for t=1:time_point_num
    v(j,t) = kinetic_param(j);

    for i=1:size(S_f_b,1)
        if S_f_b(i,j) < 0
            v(j,t)=v(j,t)*conc(i,t);
        end
    end
 end
end

%---------- Compute Vnet ----------------
Vnet = zeros(size(expanded_model.unexp_rxn_info,1),time_point_num);

for j=1:size(expanded_model.unexp_rxn_info,1)
  for t=1:time_point_num
     % If this is not an export reaction which is irreversible
     if expanded_model.unexp_rxn_info{j,6} ~= 3 
        % index of the first elementary rxn corresponding to j
        elem_index_first = expanded_model.unexp_rxn_info{j,2};

        % Corresponding indices of v_f and v_b. The first element is index of the forward
        % rxn whereas the second element is the index of the backward rxn
        v_indices=find(expanded_model.rxn_f_bInd_rxnInd==elem_index_first);

        % Vnet = v_f - v_b
        Vnet(j,t) = v(v_indices(1),t) - v(v_indices(2),t);

        % Check if v_f - v_b are equal for elemetnary steps
        for k=expanded_model.unexp_rxn_info{j,2}:expanded_model.unexp_rxn_info{j,3}
           v_k_indices=find(expanded_model.rxn_f_bInd_rxnInd==k);
           if abs(v(v_k_indices(1),t) - v(v_k_indices(2),t) - Vnet(j,t)) > 0.1
%                 fprintf('(v_f - v_b) not equal for all elementary steps of reaction %s at time %d\n',expanded_model.unexp_rxn_info{j,1},t);
             end
        end

     % if this an export rxn which is irreversilbe
     else
        % Vnet = v
        elem_index = expanded_model.unexp_rxn_info{j,2};
        v_index=find(expanded_model.rxn_f_bInd_rxnInd==elem_index);
        Vnet(j,t) = v(v_index,t);
     end
  end   % end for t
end

