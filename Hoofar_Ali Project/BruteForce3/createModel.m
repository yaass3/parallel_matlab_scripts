function [kinetic_param,elem_rate,sampled_R,sampled_R_vec,sampled_e,sampled_e_vec,MyRev] = createModel(expanded_model)

% This function creates a model for the ensumble
%
% INPUTS:
% -------
%     expanded_model: Expanded model
%   
% OUTPUTS:
% --------
%     kinetic_param: Kinetic parameters
%         elem_rate: Rate of elementary reactions
%         sampled_R: A matrix containing the sampled reversibilities for each reaction 
%                    Columns correspond to reactions in the unexpanded model plus reactions
%                    showing regulatory constraints in the network 
%                    Rows correspond to reactions in the expanded model where  
%                    reversible reactions have not been decomposed into forward and backwqrd
%     sampled_R_vec: A vector whose size is equal to the size of expanded_model.rxn
%                    and contains the reaction reversibilties
%         sampled_e: A matrix containing the enzyme fractions for each enzyme
%                    Columns correspond to enzymes. Non-zero elements in each column 
%                    correspond to enzyme fractions related to that enzyme
%                    Rows correspond metabolites in the expanded network 
%     sampled_e_vec: A vector whose size is equal to the number of metabolites in 
%                    expanded network and stores the enzyme fractions 
%
% Ali R. Zomorrodi - Costas Maranas Lab, April 2012
%

% Sample enzyme fractions (e)
[sampled_e,sampled_e_vec] = sample_e(expanded_model);

% sample reaction reversibilities
[sampled_R,sampled_R_vec,MyRev] = sample_R(expanded_model);

% Compute rate of elementary reactions
elem_rate = compute_elemRates(expanded_model,sampled_R_vec);

% Compute the kinetic parameters
kinetic_param = compute_kineticParam(expanded_model,sampled_e_vec,elem_rate);      

