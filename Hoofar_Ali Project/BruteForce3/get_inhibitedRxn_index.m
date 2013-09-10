function  inhibiied_rxn_index = get_inhibitedRxn_index(expanded_model,unexp_rxn_index,i);

% Function to find out the index of inhibited rxn 
% The index oc inhibited rxn is in expanded_model.S_f_b
%
%     unexp_rxn_index: Index of the regulatory rxn (in the unexpanded network)
%                   i: Index of the elementary rxn corresponding to unexp_rxn_index
% inhibiied_rxn_index: Index of the inhibited reaction

%---- Competetive inhibition -----
if expanded_model.unexp_rxn_info{unexp_rxn_index,6} == -1

   % find which reactant of this elementary rxn is an enzyme
   reactant_indices = find(expanded_model.S(:,i) < 0);

   % Index of enzyme
   enz_ind = reactant_indices(find(reactant_indices >= expanded_model.metabIndexInfo{2,2}));

   % Determine in which rxns this enzyme participates as a reactant
   enz_react_rxn_indices=find(expanded_model.S(enz_ind,:) < 0);

   % Find in which ones of these rxns are not regulatory
   % We consider only one of these non-regulatory reactions though
   for k=1:length(enz_react_rxn_indices)
      if expanded_model.unexp_rxn_info{get_Vnet_index(expanded_model,enz_react_rxn_indices(k)),6} > 0
         rxn_ind_nonReg(k) = enz_react_rxn_indices(k);
      end
   end

   % Find the index of this rxn in rxn_f_b 
   % Consider only one reaction in rxn_ind_nonReg e.g., the first one
   v_f_b_indices = find(expanded_model.rxn_f_bInd_rxnInd == rxn_ind_nonReg(1));
   inhibiied_rxn_index = v_f_b_indices(1);  % The first one if the forward

%---- Uncompetetive inhibition ----
elseif expanded_model.unexp_rxn_info{unexp_rxn_index,6} == -2

   % find which reactant of this elementary rxn is an enzyme complex
   reactant_indices = find(expanded_model.S(:,i) < 0);

   % Index of enzyme complex
   enzComp_ind = reactant_indices(find(reactant_indices >= expanded_model.metabIndexInfo{3,2}));

   % Determine in which rxns this enzyme complex participates as a reactant
   enzComp_react_rxn_indices=find(expanded_model.S(enzComp_ind,:) < 0);

   % Find in which ones of these rxns are not regulatory
   % We consider only one of these non-regulatory reactions though
   for k=1:length(enzComp_react_rxn_indices)
      if expanded_model.unexp_rxn_info{get_Vnet_index(expanded_model,enzComp_react_rxn_indices(k)),6} > 0
         rxn_ind_nonReg(k) = enzComp_react_rxn_indices(k);
      end
   end

   % Find the index of this rxn in rxn_f_b 
   % Consider only one reaction in rxn_ind_nonReg e.g., the first one
   v_f_b_indices = find(expanded_model.rxn_f_bInd_rxnInd == rxn_ind_nonReg(1));
   inhibiied_rxn_index = v_f_b_indices(1);  % The first one if the forward

%----- Mixed inhibition ----
elseif expanded_model.unexp_rxn_info{unexp_rxn_index,6} == -3
   % Competetive inhibition part
   if i == expanded_model.unexp_rxn_info{unexp_rxn_index,2}

     % find which reactant of this elementary rxn is an enzyme
     reactant_indices = find(expanded_model.S(:,i) < 0);

     % Index of enzyme
     enz_ind = reactant_indices(find(reactant_indices >= expanded_model.metabIndexInfo{2,2}));

     % Determine in which rxns this enzyme participates as a reactant
     enz_react_rxn_indices=find(expanded_model.S(enz_ind,:) < 0);

     % Find in which ones of these rxns are not regulatory
     % We consider only one of these non-regulatory reactions though
     for k=1:length(enz_react_rxn_indices)
        if expanded_model.unexp_rxn_info{get_Vnet_index(expanded_model,enz_react_rxn_indices(k)),6} > 0
           rxn_ind_nonReg(k) = enz_react_rxn_indices(k);
        end
     end

     % Find the index of this rxn in rxn_f_b 
     % Consider only one reaction in rxn_ind_nonReg e.g., the first one
     v_f_b_indices = find(expanded_model.rxn_f_bInd_rxnInd == rxn_ind_nonReg(1));
     inhibiied_rxn_index = v_f_b_indices(1);  % The first one if the forward

   
   % Uncompetetive inhibition part
   elseif i == expanded_model.unexp_rxn_info{unexp_rxn_index,3}

     % find which reactant of this elementary rxn is an enzyme complex
     reactant_indices = find(expanded_model.S(:,i) < 0);

     % Index of enzyme complex
     enzComp_ind = reactant_indices(find(reactant_indices >= expanded_model.metabIndexInfo{3,2}));

     % Determine in which rxns this enzyme complex participates as a reactant
     enzComp_react_rxn_indices=find(expanded_model.S(enzComp_ind,:) < 0);

     % Find in which ones of these rxns are not regulatory
     % We consider only one of these non-regulatory reactions though
     for k=1:length(enzComp_react_rxn_indices)
        if expanded_model.unexp_rxn_info{get_Vnet_index(expanded_model,enzComp_react_rxn_indices(k)),6} > 0
           rxn_ind_nonReg(k) = enzComp_react_rxn_indices(k);
        end
     end

     % Find the index of this rxn in rxn_f_b 
     % Consider only one reaction in rxn_ind_nonReg e.g., the first one
     v_f_b_indices = find(expanded_model.rxn_f_bInd_rxnInd == rxn_ind_nonReg(1));
     inhibiied_rxn_index = v_f_b_indices(1);  % The first one if the forward

   end

% Activation
elseif expanded_model.unexp_rxn_info{unexp_rxn_index,6} == -4
   fprintf('\nPlease complete part of the code ''compute_elemRates.m'' for activation !!\n');
end


