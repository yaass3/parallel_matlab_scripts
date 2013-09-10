function  elem_rate = compute_elemRates(expanded_model,sampled_R_vec)

% This function computes the rate of elementray reactions using
% eqns 12 and 13 of Trna's paper (PMID: 18820235)
%
% INPUTS:
% -------
%  expanded_model: Expanded metabolic network
%   sampled_R_vec: A vector whose size is equal to the size of expanded_model.rxn
%                  and contains the reaction reversibilties
%
% OUTPUTS:
% --------
%      elem_rate:  Rate of elementary reactions
%


elem_rate=zeros(size(expanded_model.rxn_f_b));

counter = 0;

for i=1:length(expanded_model.rxn)

   unexp_rxn_index = get_Vnet_index(expanded_model,i);

   Vnet = expanded_model.unexp_rxn_info{unexp_rxn_index,7};

   % if reversible rxn   
   if expanded_model.rev(i) == 1

      if sign(Vnet) ~= 0
        % v_forward
        counter = counter + 1;
        elem_rate(counter) = Vnet/(1-[sampled_R_vec(i)^sign(Vnet)]);

	% v_backward
        counter = counter + 1;
        elem_rate(counter) = Vnet*[sampled_R_vec(i)^sign(Vnet)]/(1-[sampled_R_vec(i)^sign(Vnet)]);

      % If Vnet is zero and this is a regulation reaction
      elseif [sign(Vnet) == 0] && [expanded_model.unexp_rxn_info{unexp_rxn_index,6} < 0]
        inhibiied_rxn_index = get_inhibitedRxn_index(expanded_model,unexp_rxn_index,i);

        % Here we assume that the flux of inhibitory reaction is a fraction of the flux of
        % inhibited reaction. We use this value instead of Vnet
        % Generate a random number between 0.01 and 0.9
        fraction = 0.01 + (0.9-0.01)*rand;

        % v_forward
        counter = counter + 1;
        elem_rate(counter) = fraction*elem_rate(inhibiied_rxn_index);

        % v_backward
        counter = counter + 1;
        elem_rate(counter) = fraction*elem_rate(inhibiied_rxn_index);

      % If Vnet is zero and this is not a regulation reaction
      elseif [sign(Vnet) == 0] && [expanded_model.unexp_rxn_info{unexp_rxn_index,6} > 0]
          fprintf('\nError! Vnet for %s is zero\n',expanded_model.unexp_rxn_info{unexp_rxn_index,1});
      end
   elseif expanded_model.rev(i) == 0
        counter = counter + 1;
        elem_rate(counter) = Vnet;
   end
end
