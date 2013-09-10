function [sampled_R,sampled_R_vec,MyRev] = sample_R(expanded_model)
%
% This function performs the sampling for enzyme fractions (e) and reversabilities (R)
%
% INPUTS:
% --------
%          model: A matlab structure containing the expanded model
%
% OUTPUTS:
% -------
%     sampled_R: A matrix containing the sampled reversibilities for each reaction 
%                Columns correspond to reactions in the unexpanded model plus reactions
%                showing regulatory constraints in the network 
%                Rows correspond to reactions in the expanded model where  
%                reversible reactions have not been decomposed into forward and backwqrd
% sampled_R_vec: A vector whose size is equal to the size of expanded_model.rxn
%                and contains the reaction reversibilties
%                
% Copyright: Ali R. Zomorrodi, Costas Maranas Lab, April 2012 
%

% Initialize sampled_e
sampled_R=zeros(length(expanded_model.rxn),size(expanded_model.unexp_rxn_info,1));

%------------------ Sample reversibilities ----------------------
% The product of R (universal gas constant) and T (Temprature)
RT=1.987*298/1000;

%++++++++++++++++++
GibbsBound=[];
MyRev=[];
%++++++++++++++++++

for i=1:size(expanded_model.unexp_rxn_info,1)  % Loop over rxns in the unexpanded model
  % Not an export rxn which is irreversible and do not consider regulaiton reactions
  if (expanded_model.unexp_rxn_info{i,6} ~= 3) && (expanded_model.unexp_rxn_info{i,6} > 0)
    dG_RT_min = min(cell2mat(expanded_model.unexp_rxn_info(i,4))/RT,cell2mat(expanded_model.unexp_rxn_info(i,5))/RT);
    dG_RT_max = max(cell2mat(expanded_model.unexp_rxn_info(i,4))/RT,cell2mat(expanded_model.unexp_rxn_info(i,5))/RT);

    if dG_RT_min == dG_RT_max
       dG_RT_min = dG_RT_min - 0.001*dG_RT_min;
       dG_RT_max = dG_RT_max + 0.001*dG_RT_max;
    end

    if sign(cell2mat(expanded_model.unexp_rxn_info(i,7))) == 1       % Sign(Vnet) = 1
      LB_R_product = min(exp(dG_RT_min),exp(dG_RT_max));
      UB_R_product = max(exp(dG_RT_min),exp(dG_RT_max));
    elseif sign(cell2mat(expanded_model.unexp_rxn_info(i,7))) == -1  % Sign(Vnet) = -1
      LB_R_product = min(exp(-dG_RT_min),exp(-dG_RT_max));
      UB_R_product = max(exp(-dG_RT_min),exp(-dG_RT_max));
    elseif sign(cell2mat(expanded_model.unexp_rxn_info(i,7))) == 0   % Sign(Vnet) = 0
      LB_R_product = 0.000001; 
      UB_R_product = 0.001; 
    end

 
    % Randomly sample between LB_R_product and UB_R_product
    % Note that R1*R2*R3 = G   and LB_R_product <=  G <= UB_R_product
    % Since R1, R2 and R3 are between zero and one G should be also between zero and one
    if LB_R_product < 0
        Lower_G_value = 0;
    else
        Lower_G_value = max(0,LB_R_product);
    end

    if UB_R_product > 1
        upper_G_value = 1;
    else
        upper_G_value = min(1,UB_R_product);
    end

    % randomize seed 
    rand('twister',sum(100*clock));

    % Generate a random number between Lower_G_value and upper_G_value
    G_sample = Lower_G_value + (upper_G_value-Lower_G_value).*rand;

    % Determine the total number of elementary steps for this reaction 
    elem_step_num = cell2mat(expanded_model.unexp_rxn_info(i,3))-cell2mat(expanded_model.unexp_rxn_info(i,2))+1;

    % Randomly sample R_prime between 0 and 1
    R_prime_sample=rand(1,elem_step_num);
    R_prime_sample = R_prime_sample./repmat(sum(R_prime_sample,2),1,elem_step_num);

    % Compute the power = log(G_sample)/log(prod(R_prime_sample))
    power_value = log(G_sample)./log(prod(R_prime_sample,2));

    % Compute R = R_prime^power_value
    R_sample = R_prime_sample.^repmat(power_value,1,elem_step_num);

    %---- Solve an optimization problem to get the value of R ---
    % initial condition
    R0 = R_sample;
  
    % Lower and upper bounds
    LB = 0.01*ones(elem_step_num,1);
    UB = 0.98*ones(elem_step_num,1);

    options = optimset('Algorithm','interior-point','display','off');

    [Ropt,fopt,exit_tag] = fmincon(@(R)objfun(R),R0,[],[],[],[],LB,UB,@(R)nonlcon(R,LB_R_product,UB_R_product),options);

    % Indices of corresponding elementary rxns
    indices=cell2mat(expanded_model.unexp_rxn_info(i,2)):cell2mat(expanded_model.unexp_rxn_info(i,3));

    
    %++++++++++++++++++++++++++
    
     GibbsBound(i,1)=log(Lower_G_value);
     GibbsBound(i,2)=log(upper_G_value);
    MyRev(i)=log(prod(Ropt));
    %++++++++++++++++++++++++++

    if exit_tag == 1
     sampled_R(indices,i) = Ropt;
    else
      fprintf('\nNo optimal solution found\n');
    end
  end
end

% Store all sampled reversabilities in a vector whose size is equal to the size of 
% exapnded_model.rxn
sampled_R_vec = zeros(length(expanded_model.rxn),1); 

for j=1:size(sampled_R,2)
   indices = find(sampled_R(:,j) ~= 0); 
   sampled_R_vec(indices) = sampled_R(indices,j);
end
% save('MyGibbs.mat','GibbsBound','MyRev');% Ali

%------- Function defining the objective function -----------
function f = objfun(R)
f=0;

%------- Function defining the constraints -----------
function [c, ceq] = nonlcon(R,LB_R_product,UB_R_product)
c(1) = -prod(R) + LB_R_product;     
c(2) = prod(R) - UB_R_product;
ceq = [];


