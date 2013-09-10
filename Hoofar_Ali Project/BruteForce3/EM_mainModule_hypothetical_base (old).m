%------------------------------------------------------------
% This code perfomrs the EM procedures
% 
% Ali R. Zomorrodi, March 2012, Costas Marans Group, 
%
%------------------------------------------------------------

clear
clc

load GlycVI_myPTS_ru5p_noHXK.mat

% Reformat the model provided by Jimmy to be consistent with the codes written by Ali
model=reformatModel(Net,[1 32],[29 30 31],'GlycVI_myPTS_ru5p_noHXK');%model=reformatModel(Net,1,'GlycVI_myPTS_ru5p_noHXK');%old version

% Expand the network
expanded_model=expandNet(model);

% Create the Jacobian matrix (Run this only once for each model)
% Use the following lines only ONCE for each input model
%fprintf('\nCreating the Jacobian matrix ...\n');
%[jacobFun_sym,f_sym,v_sym,y_vec_sym] = create_jacobian(expanded_model);


% Time interval for all ODE simulations
t_interval = [0:100];

% Threshold for accepting a model. For the hypothetical base model use a more
% stringent threshold of Vnet_thr/10
Vnet_thr = 0.15;

%---------------------------------------------------------------
%---------- Create a base hypothetical model --------------------
%---------------------------------------------------------------
fprintf('\n----- Creating a base hypothetical model -----\n');

done = 0;

% Keep creating a hypthetical model base until the output is satisfactory
while done == 0

  [kinetic_param_base,~,~,~,~,sampled_e_vec_base] = createModel(expanded_model);

  %--- Solve ODEs and system of equations for the wild type ---
  % This is to check if all models reach the same steady state
  % First create the initial condition for concentrations
  y0=sampled_e_vec_base;

  % For cofactors the initial concentraiton should be one
  y0(expanded_model.cof_metab_index) = 1;
  y0(expanded_model.cof_metab_index+1) = 1;

  % Initial normalized concentrations for other metabolites in the unexpanded model
  y0(1:expanded_model.cof_metab_index(1) -1) = rand;%Ali %rand(1:expanded_model.cof_metab_index(1) -1,1);
 
  % Solve the system of ODEs
  [T_base,conc_base,v_base,Vnet_base] = solve_ode(expanded_model,t_interval,y0,kinetic_param_base);

  % Solve system of equations
  SolveSystemOfEqn = 0;
 if SolveSystemOfEqn == 1
     fprintf('    Solving the system of equations ...\n');
    perturbation = expanded_model.Perturbations(:,1);
    [conc_ss_base,v_ss_base,Vnet_ss_base] = solve_ss(expanded_model,y0,kinetic_param_base,perturbation);
  end

  % We use one of the export fluxes (the first of them) as the basis to accept a model
  % Index of the export fluxes
  export_indices = find(cell2mat(expanded_model.unexp_rxn_info(:,6)) == 3);  

  if  abs(Vnet_base(export_indices(1),end)  - expanded_model.unexp_rxn_info{export_indices(1),7}) < (Vnet_thr/10)*expanded_model.unexp_rxn_info{export_indices(1),7}
      done = 1;
      fprintf('A hypothetical base model was created\n\n');
  else
      fprintf('The attempt for making a base hypothetical model was not successful! Trying again ...\n');
  end
 
end

%------------------------------------------------------------------------
%----------------- Creating the ensemble of models ----------------------
%------------------------------------------------------------------------

fprintf('\n----- Creating the ensemble of models -----\n');

% Specify the number of models in the ensemble
ensemble_size = 1000;

% Number of acceptable models
acceptModel_counter=0;

% Maximum allowable number of tries for populating the ensemble
allowable_try_num = 10;

% Number of tries for populating the ensemble
try_num=1;

% The index of export reactions to be used for testing if a model is acceptable 
export_indices = find(cell2mat(expanded_model.unexp_rxn_info(:,6)) == 3);

done = 0;

while done == 0

  clear conc_ode_ens_try Vnet_ode_ens_try conc_ode_ss_ens_try Vnet_ode_ss_ens_try sampled_e_vec_ens_try kinetic_param_ens_try
  
  for  n=1:(ensemble_size - acceptModel_counter)
     fprintf('    Creating model %i ...\n',n); 
     % Create the model
     [kinetic_param,~,~,~,~,sampled_e_vec] = createModel(expanded_model);
  
     % Store the kinetic parameters and sampled enzyme fractions for each model
     sampled_e_vec_ens_try(:,n)=sampled_e_vec; 
     kinetic_param_ens_try(:,n)=kinetic_param;
  
     %---- Solve ODEs and system of equations for the wild type -----
     % This is to check if all models reach the same steady state
     % First create the initial condition for concentrations
     y0=sampled_e_vec;
  
     % For cofactors the initial concentraiton should be one
     y0(expanded_model.cof_metab_index) = 1;
     y0(expanded_model.cof_metab_index+1) = 1;
  
     % Initial normalized concentrations for other metabolites in the unexpanded model
     y0(1:expanded_model.cof_metab_index(1) -1) = rand;
   
     % Solve the system of ODEs
     [T,conc,v,Vnet] = solve_ode(expanded_model,t_interval,y0,kinetic_param);
  
     % Store results for all models in the ensemble in variables
     conc_ode_ens_try(:,:,n) = conc;           % concentrations
     Vnet_ode_ens_try(:,:,n) = Vnet;           % Vnet
     T_ens(:,n) = T;                           % Time points
     conc_ode_ss_ens_try(:,n) = conc(:,end);   % Steady state concentrations 
     Vnet_ode_ss_ens_try(:,n) = Vnet(:,end);   % Steady state net fluxes
  
     % Solve system of equations
     SolveSystemOfEqn = 0;
     if SolveSystemOfEqn == 1
       fprintf('    Solving the system of equations ...\n');
       perturbation = expanded_model.Perturbations;
       [conc_ss,v_ss,Vnet_ss] = solve_ss(expanded_model,y0,kinetic_param,perturbation);
     end
  
  end   % end of for

   %--- Test if the created models are acceptable ---
  for  n=1:(ensemble_size - acceptModel_counter)
     if  abs(Vnet_ode_ens_try(export_indices(1),end,n)  - expanded_model.unexp_rxn_info{export_indices(1),7}) < Vnet_thr*expanded_model.unexp_rxn_info{export_indices(1),7}
       acceptModel_counter = acceptModel_counter + 1; 

       % Store results for accepted models in the ensemble in variables
       sampled_e_vec_ens(:,acceptModel_counter) = sampled_e_vec_ens_try(:,n);
       kinetic_param_ens(:,acceptModel_counter) = kinetic_param_ens_try(:,n);
       conc_ode_ens(:,:,acceptModel_counter) = conc_ode_ens_try(:,:,n);    
       Vnet_ode_ens(:,:,acceptModel_counter) = Vnet_ode_ens_try(:,:,n);   
     end
  end

  if acceptModel_counter == ensemble_size
      done = 1;
      fprintf('\nThe ensemble of models was successfully created!\n');
  else
      fprintf('\n  In try #i only %i were acceptable. Trying again to create the rest ...\n',try_num,acceptModel_counter);
      try_num = try_num + 1;
   
      if try_num > allowable_try_num 
        done = 1;
        fprintf('The maximum allowable number of tries was excieeded. The number of accepted models in the ensemble is less than %i. Try increasing the value of variable ''allowable_try_num''\n',ensemble_size)
      end

  end

end    % End of while

%--- Plot the results for the wild-type ----
doPlots = 0;

if doPlots == 1
  fprintf('\n----- Plot the results for the wild type -----\n')
  T=T_ens(:,1);
  clear T_ens;
  for   n=1:ensemble_size
    figure(1);
    plot(T,conc_ode_ens(20,:,n));
    hold on;

    figure(2)
    plot(T,Vnet_ode_ens(10,:,n));
    hold on;
  end
end


% Save the perturbation results in a mat file
save(strcat('ensemble_',expanded_model.model_name,'.mat'),'expanded_model','kinetic_param_ens','sampled_e_vec_ens','kinetic_param_base','sampled_e_vec_base'); 

%-------------------------------------------------------------
%---------------- Perform the perturbations -------------------
%--------------------------------------------------------------
fprintf('\n----- Perform the perturbations -----\n');



% Save the name of perturbed enzymes in a variable
perturb_enz_name = cell(1,1); 

perturb_counter=0;

%-- Index of enzymes that shuuld not be perturbed --
no_perturb_index = [expanded_model.cof_rxn_index find(cell2mat(expanded_model.unexp_rxn_info(:,6)) > 1)'];

% Max index of internal rxns
internal_rxns_maxIndex = max(find(cell2mat(expanded_model.unexp_rxn_info(:,6)) > 1));

for i=1:size(expanded_model.Perturbations,1)
   if i > internal_rxns_maxIndex
      no_perturb_index = [no_perturb_index i]; 
   end
end

% Index of the fluxes for which experimental data is available and are used to 
% screen the models in each perturbation
screen_rxn_indices = export_indices;

%-- Add perturbations to the model --
for n_perturbType = 1:2    % Perturbation type: 1 = knockout , 2 = overexpression
% for n_enz = 1:length(expanded_model.Enzyme)
 for n_enz = 4:4

  % Perturb only enzymes not in no_perturb_index 
  if isempty(find(no_perturb_index==n_enz)) 

   Ldr=0.00001;
   Loe=2;

   expanded_model.Perturbations(:,2)=ones(size(expanded_model.Perturbations,1),1);
   if n_perturbType == 1
     expanded_model.Perturbations(n_enz,2)=Ldr;
   elseif n_perturbType == 2
     expanded_model.Perturbations(n_enz,2)=Loe;
   end

   % Create the perturbation names
   expanded_model.perturb_names = createPerturbName(expanded_model);
   
   % If the user has provided any perturbation vector
   if size(expanded_model.Perturbations,2) > 1
   
    current_ensemble_size = ensemble_size; 
   
    % Kinetic parameters for the models in the current ensemble
    kinetic_param_current_ens = kinetic_param_ens;
   
    % Initial sampled enzyme fractions for the models in the current ensemble
    sampled_e_vec_current_ens = sampled_e_vec_ens;
   
   
    % Loop over perturbations
    for np = 2:size(expanded_model.Perturbations,2)
       fprintf('\n*****\n')
   
       %--- Perform the perturbations for the base model ---
       fprintf('%s for the base model ...\n\n',expanded_model.perturb_names{np});

       [~,conc_perturb_base,v_perturb_base,Vnet_perturb_base] = perturbModel(expanded_model,kinetic_param_base,t_interval,sampled_e_vec_base,np);
   
       %--- Perform the perturbations for the models in the ensemble ---
       fprintf('%s for the models in the ensemble ...\n',expanded_model.perturb_names{np});
       [metrics_perturb,Vnet_perturb_ode_ens,conc_perturb_ode_ens] = perturbEnsemble(expanded_model,kinetic_param_current_ens,sampled_e_vec_current_ens,screen_rxn_indices,t_interval);
   
       %--- Test which one of the perturbed models are acceptable ---
       acceptModel_counter = 0;
       
       kinetic_param_current_ens_old = kinetic_param_current_ens;
       sampled_e_vec_current_ens_old =  sampled_e_vec_current_ens;
       clear kinetic_param_current_ens sampled_e_vec_current_ens
    
       for  n=1:current_ensemble_size
          if  length(find(abs(Vnet_perturb_base(screen_rxn_indices,end) - Vnet_perturb_ode_ens(screen_rxn_indices,end,n)) < 0.01)) == length(screen_rxn_indices) 
            acceptModel_counter = acceptModel_counter + 1; 
   
            kinetic_param_current_ens(:,acceptModel_counter) = kinetic_param_current_ens_old(:,n);
            sampled_e_vec_current_ens(:,acceptModel_counter) = sampled_e_vec_current_ens_old(:,n);
          end
       end

       perturb_counter = perturb_counter + 1;

       Vnet_sd_ave = metrics_perturb(1,:);
       Vnet_sd_ave_screen = metrics_perturb(2,:);
       
       fprintf('# accepted = %i , # rejected = %i\n',acceptModel_counter,current_ensemble_size - acceptModel_counter);
       fprintf('SD average = %.4f  ,   SD screen = %.4f\n',Vnet_sd_ave,Vnet_sd_ave_screen);
   
       current_ensemble_size = acceptModel_counter;
   
    end  % end of for np
   end  % end of if size(expanded_model.Perturbations,2) > 1
   
   
  end    % end of if isempty(find(Net.Vin_index == n_enz))
 end % end for n_enz
end % end for n_perturbType 


% Save the perturbation results in a mat file
save(strcat('perturb_results_',expanded_model.model_name,'.mat'),'expanded_model','perturb_enz_name','metrics_perturb','kinetic_param_ens','sampled_e_vec_ens','kinetic_param_base','sampled_e_vec_base','no_perturb_index');

