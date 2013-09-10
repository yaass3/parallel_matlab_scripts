function [sampled_e,sampled_e_vec] = sample_e(expanded_model)
%
% This function performs the sampling for enzyme fractions (e) and reversabilities (R)
%
% INPUTS:
% --------
%          model: A matlab structure containing the expanded model
%
% OUTPUTS:
% -------
%     sampled_e: A matrix containing the enzyme fractions for each enzyme
%                Columns correspond to enzymes. Non-zero elements in each column 
%                correspond to enzyme fractions related to that enzyme
%                Rows correspond metabolites in the expanded network 
% sampled_e_vec: A vector whose size is equal to the number of metabolites in 
%                expanded network and stores the enzyme fractions 
%                
% Copyright: Ali R. Zomorrodi, Costas Maranas Lab, April 2012 
%

% Initialize sampled_e
sampled_e=zeros(size(expanded_model.enz_enzComplex));

%------------------ Sample enzyme fractions ----------------------
% Loop over each enzyme
counter=0;
for j=expanded_model.metabIndexInfo{2,2}:expanded_model.metabIndexInfo{2,3}

  counter=counter+1;

  % randomize seed 
  rand('twister',sum(100*clock));

  % Total number of enzyme fractions associated with each enzyme
  % One in the RHS is the free enzyme and second term determines the number of enzyme 
  % complexes associated with that enzyme
  enz_frac_num = 1 + sum(expanded_model.enz_enzComplex(:,counter)); 

  % Perform the sampling
  e_sample_pool=rand(10000,enz_frac_num);

  % Normalize
  e_sample_pool = e_sample_pool./repmat(sum(e_sample_pool,2),1,enz_frac_num);

  % Draw one sample from the sample pool randomly
  % The index of enzyme and complxes associated with enzyme j
  indices = [j;find(expanded_model.enz_enzComplex(:,counter)~=0)];

  isEqualToOne = 0;
  while isEqualToOne == 0
      random_index=randi(10000,1);
      if sum(e_sample_pool(random_index,:)) == 1
          isEqualToOne = 1;
      end
  end
  sampled_e(indices,counter) = e_sample_pool(random_index,:);
end


% Store all enzyme fractions in a vector whose size is equal to the number of 
% metabolites in the expanded model (expanded_model.metab)
sampled_e_vec = zeros(size(expanded_model.metab));

for i=1:size(sampled_e,2)
    % Index of non-zero elements
    non_zero_indices = find(sampled_e(:,i) ~= 0);
    sampled_e_vec(non_zero_indices) = sampled_e(non_zero_indices,i); 
end

