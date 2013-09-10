function [metrics_perturb,Vnet_perturb_ode_ens,conc_perturb_ode_ens] = perturbEnsemble(expanded_model,kinetic_param_ens,sampled_e_vec_ens,screen_rxn_indices,t_interval)

% This function computes the primal problem, i.e., given a set of binary variables
% solves the nonlinear problem (to compute variance)
%
% INPUTS:
% -------
%      expanded_model: Expanded metabolic model
%  kinetic_param_ens:  kinetic_param_ens
%  sampled_e_vec_ens:  Sample e vector   
% screen_rxn_indices:  Rxn indices on which the variance is calculated
%
% OUTPUTS:
% -------
%     metrics_pertub:   A matrix containing the perturbation results
%                       The columns of the matrix are equal to the number of sets of perturbaitons
%                       (i.e., the number of columns of expanded_model.Perturbations 
%                       minus the first column)
%                       1st row: Average standard deviation over all fluxes 
%                       2nd row: Average standard deviation over the fluxes
%                       specified by screen_rxn_indices 
% Vnet_perturb_ode_ens: Net fluxes after perturbation
% conc_perturb_ode_ens: Concentrations after pertrubations     
%
% Ali R. Zomorrodi- Costas Maranas Lab, May 2012
%


% 1st row: Average standard deviation over all fluxes 
% 2nd row: Average standard deviation over the fluxes used as the screenning basis
metrics_perturb = zeros(2+length(screen_rxn_indices),1);

% Size of the ensemble
ensemble_size = size(kinetic_param_ens,2);

% Index of rxns that are not regulatory
nonReg_rxn_index = find(cell2mat(expanded_model.unexp_rxn_info(:,6)) > 0);

% Flux of the wildtype strain at steady state
wild_ss = cell2mat(expanded_model.unexp_rxn_info(nonReg_rxn_index,7));

perturb_counter=0;

% Loop over perturbations
for np = 2:size(expanded_model.Perturbations,2)

   %--- Perform the perturbations for the models in the ensemble ---

   % Loop over the models in the ensemble
   for n = 1:ensemble_size

        % Perform the perturbation
        [T,conc_perturb,v_perturb,Vnet_perturb] = perturbModel(expanded_model,kinetic_param_ens(:,n),t_interval,sampled_e_vec_ens(:,n),np);

        % Store results for all models in the ensemble in variables
        conc_perturb_ode_ens(:,:,n) = conc_perturb;    
        Vnet_perturb_ode_ens(:,:,n) = Vnet_perturb;   
        T_ens(:,n) = T;                              

   end   % end of parfor n

   Vnet_perturb_ss_ode_ens(:,:) = Vnet_perturb_ode_ens(:,end,:);

   % Mean of each reaction over all models = sum(n,Vnet(j,n)
   Vnet_mean_ens = mean(Vnet_perturb_ss_ode_ens,2);

   % Standard deviation for each reaction over all models in the ensemble
   % sqrt(sum(n,Vnet(j,n) - Vnet_mean)/(N_ensemble - 1))
   %Vnet_sd_ens=sqrt(sum([Vnet_perturb_ss_ode_ens - repmat(Vnet_mean_ens,1,ensemble_size)].^2,2)/(ensemble_size - 1));
   Vnet_sd_ens=sum([Vnet_perturb_ss_ode_ens - repmat(Vnet_mean_ens,1,ensemble_size)].^2,2)/(ensemble_size - 1);

   % Sum of the standard deviation for all reactions 
   % Since the numbers are usually small multiply them by 100 to make the comparison 
   % easier
   %Vnet_sd_ave = 100*sum(Vnet_sd_ens)/length(Vnet_sd_ens);
   Vnet_sd_ave = sum(Vnet_sd_ens)/length(Vnet_sd_ens);

   %--- Compute the average variance only for screening rxns  ---
   % Screening reactions are those for which experimental data is avaiable and are
   % used to screen the models after each perturbation
   %Vnet_sd_ave_screen = 100*sum(Vnet_sd_ens(screen_rxn_indices))/length(screen_rxn_indices);
   UPEX=[];
   %++++++++++++++++++++
   for k=1:length(screen_rxn_indices)
       UPEX(k)=Vnet_sd_ens(screen_rxn_indices(k));
   end
   %++++++++++++++++++++
   Vnet_sd_ave_screen = sum(Vnet_sd_ens(screen_rxn_indices))/length(screen_rxn_indices);

   perturb_counter = perturb_counter + 1;
   metrics_perturb(1,perturb_counter) = Vnet_sd_ave;
   metrics_perturb(2,perturb_counter) = Vnet_sd_ave_screen;
   metrics_perturb(3:end,perturb_counter) =UPEX;

end  % end of for np

