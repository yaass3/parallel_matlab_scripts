function expanded_model = expandNetg(model);
% This function computes the expanded stoichiometric matrix
%
%  INPUT
%---------
%           metab: Name of the metabolites
%             rxn: Name of the rxns 
%               S: Stoichiomteric matrix        
%          Enzyme: Name of enzymes
%       enzymeReg: A matrix represenitng which metabolite affect the activity of which
%                  enzyme. Rows correspond to metabolites and columns correspond to enzymes
%                  If element (i,j) of this matrix is non-zeor it implies that metabolite
%                  if has an inhibiton/activation effect on enzyme j:
%                  -1 = competetive inhibution  -2 = Uncompetetive inhibition
%                  -3 = mixed inhibution        -4 = activation 
%        rxntype:  Type of rxn: 1=internal  2=uptake 3=export
%            Vss:  Steady state fluxes. The first column correspond to the wild-type.
%                  The rest of the columns correspond to perturbations 
%            GPR:  A matrix representing the gene-(protein)-reaction associatons.
%                  Rows correspond to genes and columns to reactions
%           SGFE:  Standard Gibbs free energy
%     metabRange:  Range of metabolite concentrations                  
%      cof_metab_index:  Index oc the cofators
%  perturb_names:  Name of the perturbed strain. The first element is wild showing the
%                  wild type strain (i.e., no perturbations)     
%
% OUTPUT
%-------
% The output is a MATLAB structure containing the folloiwing fields
%           metab:  Name of the metabolites in the expanded model
%             rxn:  Name of the rxns in the expanded model, where all elementtary rxns are
%                   considered as reversible
%         rxn_f_b:  Name of the rxns in the expanded model where all elementary rxns 
%                   are decomposed into two spearate forward and backward rxn
%               S:  Stoichiomteric matrix of the expanded model where the elemetnary rxns 
%                   are considered as reversible
%           S_f_b:  Stoichiomteric matrix of the expanded model where the elementary rxns 
%                   are decomposed into forward and backward 
%   enz_enzComplex: A matrix representing which enzyme complex is associated with 
%                   which enzyme. Rows correspond the metabolites of the extended model
%                   where columns correspond to enzymes. Element (i,j) of this matrix
%                   is one if an enzyme complex i is associated with enzyme j and it 
%                   is zero otherwise
% metab_enzComplex: A matrix representing which enzyme complexes are associated with 
%                   each metabolite in the unexpanded model. Rows correspond to the 
%                   metabolites of the expanded model whereas columns correspond to the
%                   metabolites of the unexpanded model. Element (i,j) of thi smatrix is
%                   one if enzyme complex i is assoicated with metabolite j in the 
%                   unexpanded model and zero otherwise
%   metabIndexInfo: Contains the start and end index of the metabolites in the unxpanded
%                   model, enzymes and enzyme complexes. Here is an example: 
%                   'Unexpanded'        [ 1]    [ 4]
%                   'Enzyme'            [ 5]    [ 9]
%                   'Enzyme complex'    [10]    [19]
%  cof_metab_index: Index of the cofators
%              rev: 1 if a reaction (in the filed rxn) is reversible and zero otherwise
%           LB_rxn: Lower bound on rxns stored in the filed rxn
%           UB_rxn: Upper bound on rxns stored in the filed rxn
%     isElementary: 1 if a reaction (stored in the field rxn) is elementary and 0 otherwise
%   unexp_rxn_info: Information about the rxns of the unexapnded model
%                   1st column: The name of rxn
%                   2nd and 3rd columns: The start and end indices of the corresponding 
%                   elementary rxns. 
%                   4th and 5th columns: The min and max of delta G of each reaction, 
%                   respectively
%                   6th column: The type of rxn: 1=internal  2=uptake 3=export
%                   -1 = competetive inhibition  -2 = Uncompetetive inhibition
%                   -3 = mixed inhibition  -4 = activation
%                   7th column: The steady-state flux of each rxn for the referecen
%                   steady-state (wild-type)
%                   8th, 9th, ...  columns: The steady-states of the perturbed strains
%rxn_f_bInd_rxnInd: Shows what the new index of each reaction in the field rxn is
%                   in rxn_f_b
%    perturb_names: Name of the perturbed strain. The first element is wild showing the
%                   wild type strain (i.e., no perturbations)     
%       model_name: Name of the model
% unexpanded_model: The whole unexpanded model
%
%
% Copyright: Ali R. Zomorrodi, April 2012
%            Chemical & Biological Systems Optimization Lab @ Penn State
%                 

% Number of reactions in the original model
rxn_num_orig = size(model.S,2);  

% Number of metabolites in the original model
metab_num_orig = size(model.S,1);

% Find out the number of enzymes in the network
enzyme_num = max(size(model.Enzyme));

% Number of the metabolites in expanded_S
% The first metab_num_orig columns of the expanded_S correspond to the metabolites present
% in the original network. The next enzyme_num column corressponds to the enzymes and the 
% following enzyme_complex_num columns correspond to the enzyme complexes
% We don't know how many enzyme complexes we have at this point so we consider only 
% the metabolites of the unexpanded model and 
expanded_S_metabs_num = metab_num_orig + enzyme_num;

%------------- Create the metabolite names for the expanded model -------------------
% Create a null cell array to store the name of metabolites in the expanded model
expanded_model.metab=cell(expanded_S_metabs_num,1);

% The following field stores the infomration about the index of metabolites
% The first column represents the type of metabolite:
% Original: The metabolite names in the original network
% Enzyme: Enzymes
% Enzyme complex: Enzyme complexes
% The second column represents the start index and the third column the end index
expanded_model.metabIndexInfo = cell(3,3);

% The first metab_num_orig elements correspond to the metabolites in the original model
for i=1:metab_num_orig
   expanded_model.metab(i) = model.metab(i);
end

expanded_model.metabIndexInfo(1,:) = {'Unexpanded' 1 metab_num_orig};

% Now add enzymes to the list of metabolites 
for i=1:enzyme_num
   expanded_model.metab(i + metab_num_orig) = model.Enzyme(i);
end

expanded_model.metabIndexInfo(2,:) = {'Enzyme' metab_num_orig+1 metab_num_orig+enzyme_num};

%----- Create a matrix representing which enzymes is related to which enzyme complexes --------
% The columns in this matrix correspond to enzymes and the rows correspond to metabolite names 
expanded_model.enz_enzComplex=zeros(expanded_S_metabs_num,enzyme_num);

%i---- Create a matrix representing which metabolites are related to which enzyme complexes ----
% Create a matrix representing which enzymes related to which enzyme complexes
expanded_model.metab_enzComplex=zeros(expanded_S_metabs_num,metab_num_orig);

%------------- Create the rxn names and related info for the expanded model -------------
% Create a null cell array to store the name of rxns in the expanded model
expanded_model.rxn=cell(1,1);

% Create a null cell array for stoichiomteric matrix
expanded_model.S=zeros(expanded_S_metabs_num,1);

% Store reversability, lowerbound and upperbound for elementary reactions
expanded_model.rev=zeros(1,1);
expanded_model.LB_rxn=zeros(1,1);
expanded_model.UB_rxn=zeros(1,1);

% Create a field for storing which reactions in the expanded model are elementary rxns
% and which one are non-enzymatic and hence non-elemnetary rxns
expanded_model.isElementary=zeros(1,1);

% Create a filed containing the name of perturbations
expanded_model.perturb_names=model.perturb_names;

% Number of perturbations (including the wild-type)
num_perturb = max(size(model.perturb_names));

% Storing the information related to the original rxns
% 1st column: The name of rxn
%
% 2nd and 3rd columns: The start and end indices of the corresponding 
% elementary rxns. 
%
% 4th and 5th columns: The min and max of delta G of each reaction, 
% respectively, which is computed as:
% delta G = delta G standard + 
%                        RT*ln(prod(concentration of metabolites^stoichiomteric coeff))
% But here we compute it as:
% delta G min = delta G standard + 
%        RT*ln((min conc of products^stoic coeff)/(max conc of reactants^stoic coeff))
% delta G max = delta G standard + 
%        RT*ln((max conc of products^stoic coeff)/(min conc of reactants^stoic coeff))
%
% 6th column: The type of rxn: 1=internal  2=uptake 3=export
%
% 7th column: The steady-state flux of each rxn for the referecen steady-state
%
% 8th, 9th, ...  columns: The steady-states of the perturbed strains
%
expanded_model.unexp_rxn_info=cell(rxn_num_orig,6+num_perturb);

% To compute min and max of deltaG we need the following
% Here we consider only the four metabolites that are present in the original model
Sreact=model.S;   % Contains only negative elements of S (i.e., reactants)
Sreact(Sreact>0)=0;      % set any positive elements to zero

Sprod=model.S;         % Contains only positive elements of S (i.e., products)
Sprod(Sprod<0)=0;      % Set any negative elements to zero

% The product of R (universal gas constant) and T (Temprature)
RT=1.987*298/1000;


%--------------------------  Expand the network ---------------- 
metab_index = metab_num_orig+enzyme_num;   % Counts the number of enzyme complexes
rxn_index = 0;     % Counts the number of rxns in the expanded netowrk
 
for j=1:rxn_num_orig
   % Original rxn name
   expanded_model.unexp_rxn_info(j,1)=model.rxn(j); 
   % min deltaG
   expanded_model.unexp_rxn_info(j,4)={model.SGFE(j) + RT*Sprod(:,j)'*log(model.metabRange(:,1)) + RT*Sreact(:,j)'*log(model.metabRange(:,2))}; 
   % max deltaG
   expanded_model.unexp_rxn_info(j,5)={model.SGFE(j) + RT*Sprod(:,j)'*log(model.metabRange(:,2)) + RT*Sreact(:,j)'*log(model.metabRange(:,1))}; 
   % rxntype 
   expanded_model.unexp_rxn_info(j,6)={model.rxntype(j)}; 
   % steady-state fluxes for the wildtype and perturbed strains
   for np = 1:num_perturb
      expanded_model.unexp_rxn_info(j,6+np)={model.Vss(j,np)}; 
   end 

  if sum(model.GPR(:,j)) > 0   % Only rxns which are coded by enzymes
      enzyme_ind = find(model.GPR(:,j) ~= 0);

      % Internal reactions
      if [model.rxntype(j)==1] 
         reactant_indices = find(model.S(:,j) < 0);
         product_indices = find(model.S(:,j) > 0);

         if max(size(enzyme_ind)) == 1

           %---- STEP 1 -------
           rxn_name_counter = 0;
           for k=1:length(reactant_indices)

             % Create rxn name j_k
             rxn_name_counter = rxn_name_counter + 1;
             rxn_index = rxn_index + 1;
             expanded_model.rxn(rxn_index) = strcat(model.rxn(j),'_',num2str(rxn_name_counter));

             % Create stoichiomteric matrix: S(A,J_1)=-1 and S(E,J_1)=-1
             expanded_model.S(reactant_indices(k),rxn_index) = -1;
             if k==1      % S(E,J_1) = -1
                expanded_model.S(metab_num_orig + enzyme_ind,rxn_index) = -1;
             else         % S(prev_enz_complex,j_k) = -1
                expanded_model.S(metab_index,rxn_index) = -1;
             end

             expanded_model.rev(rxn_index)=1;
             expanded_model.LB_rxn(rxn_index)=-1000;
             expanded_model.UB_rxn(rxn_index)=1000;
             expanded_model.isElementary(rxn_index)=1;
             if k==1
                expanded_model.unexp_rxn_info(j,2)={rxn_index};  % Start index
             end

             % New enzyme complex with reactants
             metab_index = metab_index + 1;
             complex_name=model.Enzyme(enzyme_ind);
             for m=1:k
                complex_name=strcat(complex_name,'_',model.metab(reactant_indices(m)));
             end
             expanded_model.metab(metab_index)=strcat(complex_name,'_complex');

             expanded_model.enz_enzComplex(metab_index,enzyme_ind) = 1;
             expanded_model.metab_enzComplex(metab_index,reactant_indices(k)) = 1;

             % S(Enzyme_complex,j_k)=1
             expanded_model.S(metab_index,rxn_index) = 1;

           end    % end for k

           %---- STEP 2 --------
           % Create rxn EAB...(other reactants) <--> ECD...(other products)
           rxn_name_counter = rxn_name_counter + 1;
           rxn_index = rxn_index + 1;
           expanded_model.rxn(rxn_index) = strcat(model.rxn(j),'_',num2str(rxn_name_counter));

           % S(EAB...,j_k) = -1;
           expanded_model.S(metab_index,rxn_index) = -1;   

           expanded_model.rev(rxn_index)=1;
           expanded_model.LB_rxn(rxn_index)=-1000;
           expanded_model.UB_rxn(rxn_index)=1000;
           expanded_model.isElementary(rxn_index)=1;

           % ECD...
           metab_index = metab_index + 1;
           complex_name=model.Enzyme(enzyme_ind);
           for m=1:length(product_indices)
              complex_name=strcat(complex_name,'_',model.metab(product_indices(m)));
           end
           expanded_model.metab(metab_index)=strcat(complex_name,'_complex');

           expanded_model.enz_enzComplex(metab_index,enzyme_ind) = 1;
           for m=1:length(product_indices)
             expanded_model.metab_enzComplex(metab_index,product_indices(m)) = 1;
           end
           % S(ECD...,J_k)=1
           expanded_model.S(metab_index,rxn_index) = 1;

           %----- STEP 3 -------
           % At each step one product is dissociated from the enzyme complex
           for k=1:length(product_indices)
             % Create rxn EAB...(other reactants) <--> ECD...(other products)
             rxn_name_counter = rxn_name_counter + 1;
             rxn_index = rxn_index + 1;
             expanded_model.rxn(rxn_index) = strcat(model.rxn(j),'_',num2str(rxn_name_counter));

             % S(previous_complex...,j_k) = -1;
             expanded_model.S(metab_index,rxn_index) = -1;

             % S(original metab k,j_k) = 1
             expanded_model.S(product_indices(k),rxn_index) = 1;
 
             expanded_model.rev(rxn_index)=1;
             expanded_model.LB_rxn(rxn_index)=-1000;
             expanded_model.UB_rxn(rxn_index)=1000;
             expanded_model.isElementary(rxn_index)=1;

           % EDF... (New complex which does not contain the first k reactants)
           if k < length(product_indices)
               metab_index = metab_index + 1;
               complex_name=model.Enzyme(enzyme_ind);
               for m=k+1:length(product_indices)
                  complex_name=strcat(complex_name,'_',model.metab(product_indices(m)));
               end
               expanded_model.metab(metab_index)=strcat(complex_name,'_complex');
 
               expanded_model.enz_enzComplex(metab_index,enzyme_ind) = 1;
               for m=k+1:length(product_indices)
                 expanded_model.metab_enzComplex(metab_index,product_indices(m)) = 1;
               end
               % S(EDF...,J_k)=1
               expanded_model.S(metab_index,rxn_index) = 1;
             else       % S(E,j_k) = 1
               expanded_model.S(metab_num_orig+enzyme_ind,rxn_index) = 1;
             end     % end of if k 

           end    % end for k

           % End index
           expanded_model.unexp_rxn_info(j,3)={rxn_index}; 

         end    % if max(size(enzyme_ind)) == 1

      % Exchnage uptake rxn  (A_ex) --> A    E --> E_A_ex   E_A_ex <--> E_A <--> E + A
      elseif (model.rxntype(j)==2) && [max(size(find(model.S(:,j) < 0))) == 0] && [max(size(find(model.S(:,j) > 0))) == 1] 
         if max(size(enzyme_ind)) == 1
           product_indices = find(model.S(:,j) > 0);

           % Create rxn name j_1: E <--> E_A_ex
           rxn_index = rxn_index + 1;
           expanded_model.rxn(rxn_index) = strcat(model.rxn(j),'_1');

           expanded_model.rev(rxn_index)=1;
           expanded_model.LB_rxn(rxn_index)=-1000;
           expanded_model.UB_rxn(rxn_index)=1000;
           expanded_model.isElementary(rxn_index)=1;
           expanded_model.unexp_rxn_info(j,2)={rxn_index};  % Start index

           % S(E,J_1)=-1
           expanded_model.S(metab_num_orig + enzyme_ind,rxn_index) = -1;

           % Write E_A_ex
           metab_index = metab_index + 1;
           expanded_model.metab(metab_index)=strcat(model.Enzyme(enzyme_ind),'_',model.metab(product_index),'_ext_complex');

           expanded_model.enz_enzComplex(metab_index,enzyme_ind) = 1;

           % S(E_A_ext,j_1)=1
           expanded_model.S(metab_index,rxn_index) = 1;

           % Create rxn name j_2: E_A_e <--> EA
           rxn_index = rxn_index + 1;
           expanded_model.rxn(rxn_index) = strcat(model.rxn(j),'_2');

           expanded_model.rev(rxn_index)=1;
           expanded_model.LB_rxn(rxn_index)=-1000;
           expanded_model.UB_rxn(rxn_index)=1000;
           expanded_model.isElementary(rxn_index)=1;
         
           % S(E_A_ext,J_2)=-1
           expanded_model.S(metab_index,rxn_index) = -1;
 
           % Write E_A_e
           metab_index = metab_index + 1;
           expanded_model.metab(metab_index)=strcat(model.Enzyme(enzyme_ind),'_',model.metab(product_index),'_complex');

           expanded_model.enz_enzComplex(metab_index,enzyme_ind) = 1;
           expanded_model.metab_enzComplex(metab_index,product_index) = 1;

           % S(E_A_e,J_2)=1
           expanded_model.S(metab_index,rxn_index) = 1;

           % Create rxn name j_3: E_A <--> E + A
           rxn_index = rxn_index + 1;
           expanded_model.rxn(rxn_index) = strcat(model.rxn(j),'_3');

           % S(E_A_e,J_3)=-1    S(E,J_3)=1   S(A_e,j_3)=1
           expanded_model.S(metab_index,rxn_index) = -1;
           expanded_model.S(metab_num_orig + enzyme_ind,rxn_index) = 1;
           expanded_model.S(product_index,rxn_index) = 1;

           expanded_model.rev(rxn_index)=1;
           expanded_model.LB_rxn(rxn_index)=-1000;
           expanded_model.UB_rxn(rxn_index)=1000;
           expanded_model.isElementary(rxn_index)=1;
           expanded_model.unexp_rxn_info(j,3)={rxn_index}; % end index
         end  % end of if max(size(enzyme_ind)) == 1

      % Exchnage uptake rxns  (A_ext) + B --> A
      elseif [model.rxntype(j)==2] && [max(size(find(model.S(:,j) < 0))) == 1] && [max(size(find(model.S(:,j) > 0))) == 1]
         reactant_indices = find(model.S(:,j) < 0);
         product_indices = find(model.S(:,j) > 0);

         % If catalyzed by only one enzyme
         if max(size(enzyme_ind)) == 1
           % Create rxn name j_1: E + (A_ext) <--> E_A_ext
           rxn_index = rxn_index + 1;
           expanded_model.rxn(rxn_index) = strcat(model.rxn(j),'_1');

           expanded_model.rev(rxn_index)=1;
           expanded_model.LB_rxn(rxn_index)=-1000;
           expanded_model.UB_rxn(rxn_index)=1000;
           expanded_model.isElementary(rxn_index)=1;
           expanded_model.unexp_rxn_info(j,2)={rxn_index};  % Start index

           % S(E,J_1)=-1
           expanded_model.S(metab_num_orig + enzyme_ind,rxn_index) = -1;

           % E_A_ext
           metab_index = metab_index + 1;
           expanded_model.metab(metab_index)=strcat(model.Enzyme(enzyme_ind),'_',model.metab(product_indices(1)),'_ext_complex');

           expanded_model.enz_enzComplex(metab_index,enzyme_ind) = 1;
           expanded_model.metab_enzComplex(metab_index,reactant_indices(1)) = 1;

           % S(E_A_ext,j_1)=1
           expanded_model.S(metab_index,rxn_index) = 1;

           % J_2 (E_A_ext + B <--> E_A_ext_B) and S(E_A_ext,J_2)=-1 and S(B,J_2)=-1
           rxn_index = rxn_index + 1;
           expanded_model.rxn(rxn_index) = strcat(model.rxn(j),'_2');
           expanded_model.S(metab_index,rxn_index) = -1;        % S(EA,j_2)=-1
           expanded_model.S(reactant_indices(1),rxn_index) = -1;  % S(B,j_2)=-1 

           expanded_model.rev(rxn_index)=1;
           expanded_model.LB_rxn(rxn_index)=-1000;
           expanded_model.UB_rxn(rxn_index)=1000;
           expanded_model.isElementary(rxn_index)=1;
 
           % E_A_ext_B
           metab_index = metab_index + 1;
           expanded_model.metab(metab_index)=strcat(model.Enzyme(enzyme_ind),'_',model.metab(product_indices(1)),'_ext_',model.metab(reactant_indices(1)),'_complex');

           expanded_model.enz_enzComplex(metab_index,enzyme_ind) = 1;
           expanded_model.metab_enzComplex(metab_index,reactant_indices(1)) = 1;

           % S(E_A_ext_B,j_2)=1
           expanded_model.S(metab_index,rxn_index) = 1;        % S(EAB,j_2)=1

           % Create rxns J_3 (E_A_ext_B <--> E_A_e) and S(E_A_ext_b,J_3)=-1 
           rxn_index = rxn_index + 1;
           expanded_model.rxn(rxn_index) = strcat(model.rxn(j),'_3');
           expanded_model.S(metab_index,rxn_index) = -1;        % S(E_A_ext_B,j_3)=-1

           expanded_model.rev(rxn_index)=1;
           expanded_model.LB_rxn(rxn_index)=-1000;
           expanded_model.UB_rxn(rxn_index)=1000;
           expanded_model.isElementary(rxn_index)=1;

           % E_A_e
           metab_index = metab_index + 1;
           expanded_model.metab(metab_index)=strcat(model.Enzyme(enzyme_ind),'_',model.metab(product_indices),'_complex');

           expanded_model.enz_enzComplex(metab_index,enzyme_ind) = 1;
           expanded_model.metab_enzComplex(metab_index,product_indices) = 1;

           % S(EC,j_3)=1
           expanded_model.S(metab_index,rxn_index) = 1;        % S(EC,j_3)=1

           % Create rxns J_4 (E_A_e <--> E + A_e) and S(E_A_e,J_4)=-1 and S(E,J_4)=1
           % and S(A_e,J_4)=1
           rxn_index = rxn_index + 1;
           expanded_model.rxn(rxn_index) = strcat(model.rxn(j),'_4');
           expanded_model.S(metab_index,rxn_index) = -1;                % S(EC,j_4)=-1
           expanded_model.S(metab_num_orig + enzyme_ind,rxn_index) = 1;   % S(E,j_4)=1
           expanded_model.S(product_indices,rxn_index) = 1;               % S(E,j_2)=1 

           expanded_model.rev(rxn_index)=1;
           expanded_model.LB_rxn(rxn_index)=-1000;
           expanded_model.UB_rxn(rxn_index)=1000;
           expanded_model.isElementary(rxn_index)=1;
           expanded_model.unexp_rxn_info(j,3)={rxn_index}; % end index

         end   % end of if max(size(enzyme_ind)) == 1

      % Exchange uptake rxn  (A_ext) + B --> A_e + C
      elseif [model.rxntype(j)==2] && [max(size(find(model.S(:,j) < 0))) == 1] && [max(size(find(model.S(:,j) > 0))) == 2]
         reactant_indices = find(model.S(:,j) < 0);
         product_indices = find(model.S(:,j) > 0);

         if max(size(enzyme_ind)) == 1

           % Create rxn name j_1: (A_ext) + E <--> E_A_ext
           rxn_index = rxn_index + 1;
           expanded_model.rxn(rxn_index) = strcat(model.rxn(j),'_1');

           % Create stoichiomteric matrix: S(E,J_1)=-1
           expanded_model.S(metab_num_orig + enzyme_ind,rxn_index) = -1;

           expanded_model.rev(rxn_index)=1;
           expanded_model.LB_rxn(rxn_index)=-1000;
           expanded_model.UB_rxn(rxn_index)=1000;
           expanded_model.isElementary(rxn_index)=1;
           expanded_model.unexp_rxn_info(j,2)={rxn_index};  % Start index

           % E_A_ext
           metab_index = metab_index + 1;
           expanded_model.metab(metab_index)=strcat(model.Enzyme(enzyme_ind),'_',model.metab(product_indices(1)),'_ext_complex');

           expanded_model.enz_enzComplex(metab_index,enzyme_ind) = 1;
           expanded_model.metab_enzComplex(metab_index,reactant_indices(1)) = 1;

           % S(E_A_ext,j_1)=1
           expanded_model.S(metab_index,rxn_index) = 1;

           % Create rxns J_2 (E_A_ext + B <--> E_A_ext_B) and S(E_A_ext,J_2)=-1 and S(B,J_2)=-1
           rxn_index = rxn_index + 1;
           expanded_model.rxn(rxn_index) = strcat(model.rxn(j),'_2');
           expanded_model.S(metab_index,rxn_index) = -1;          % S(E_A_ext,j_2)=-1
           expanded_model.S(reactant_indices(1),rxn_index) = -1;  % S(B,j_2)=-1 

           expanded_model.rev(rxn_index)=1;
           expanded_model.LB_rxn(rxn_index)=-1000;
           expanded_model.UB_rxn(rxn_index)=1000;
           expanded_model.isElementary(rxn_index)=1;
 
           % E_A_ext_B
           metab_index = metab_index + 1;
           expanded_model.metab(metab_index)=strcat(model.Enzyme(enzyme_ind),'_',model.metab(product_indices(1)),'_ext_',model.metab(reactant_indices(1)),'_complex');

           expanded_model.enz_enzComplex(metab_index,enzyme_ind) = 1;
           expanded_model.metab_enzComplex(metab_index,reactant_indices(1)) = 1;


           % S(E_A_ext_B,j_2)=1
           expanded_model.S(metab_index,rxn_index) = 1;      

           % Create rxns J_3 (E_A_ext_B <--> E_A_e_C) and S(E_A_ext_B,J_3)=-1 
           rxn_index = rxn_index + 1;
           expanded_model.rxn(rxn_index) = strcat(model.rxn(j),'_3');
           expanded_model.S(metab_index,rxn_index) = -1;        % S(E_A_ext_B,j_3)=-1

           expanded_model.rev(rxn_index)=1;
           expanded_model.LB_rxn(rxn_index)=-1000;
           expanded_model.UB_rxn(rxn_index)=1000;
           expanded_model.isElementary(rxn_index)=1;

           % E_A_e_C
           metab_index = metab_index + 1;
           expanded_model.metab(metab_index)=strcat(model.Enzyme(enzyme_ind),'_',model.metab(product_indices(1)),'_',model.metab(product_indices(2)),'_complex');

           expanded_model.enz_enzComplex(metab_index,enzyme_ind) = 1;
           expanded_model.metab_enzComplex(metab_index,product_indices(1)) = 1;
           expanded_model.metab_enzComplex(metab_index,product_indices(2)) = 1;

           % S(ECC,J_3)=1
           expanded_model.S(metab_index,rxn_index) = 1;

           % Create rxn name j_4: E_A_e_C <--> EC + A_e
           rxn_index = rxn_index + 1;
           expanded_model.rxn(rxn_index) = strcat(model.rxn(j),'_4');

           expanded_model.rev(rxn_index)=1;
           expanded_model.LB_rxn(rxn_index)=-1000;
           expanded_model.UB_rxn(rxn_index)=1000;
           expanded_model.isElementary(rxn_index)=1;

           % S(E_A_e_C,J_4) = -1
           expanded_model.S(metab_index,rxn_index) = -1;

           % EC
           metab_index = metab_index + 1;
           expanded_model.metab(metab_index)=strcat(model.Enzyme(enzyme_ind),'_',model.metab(product_indices(2)),'_complex');

           expanded_model.enz_enzComplex(metab_index,enzyme_ind) = 1;
           expanded_model.metab_enzComplex(metab_index,product_indices(2)) = 1;

           % S(EC,j_4)=1  and S(A_e,j_4) = 1
           expanded_model.S(metab_index,rxn_index) = 1;
           expanded_model.S(product_indices(1),rxn_index) = 1;

           % Create rxn name j_5: EC <--> E + C
           rxn_index = rxn_index + 1;
           expanded_model.rxn(rxn_index) = strcat(model.rxn(j),'_5');

           expanded_model.rev(rxn_index)=1;
           expanded_model.LB_rxn(rxn_index)=-1000;
           expanded_model.UB_rxn(rxn_index)=1000;
           expanded_model.isElementary(rxn_index)=1;
           expanded_model.unexp_rxn_info(j,3)={rxn_index}; % end index

           % S(EC,j_5)=-1  S(E,j5)=1  S(C,j5)=1 
           expanded_model.S(metab_index,rxn_index) = -1;
           expanded_model.S(metab_num_orig + enzyme_ind,rxn_index) = 1;
           expanded_model.S(product_indices(2),rxn_index) = 1;

         end    % if max(size(enzyme_ind)) == 1

      % Phosphotrasnferase system with glucose as external metabolite 
      % (glc_ext) + pep --> g6p + pyr (g6p should be the first metabolite and pyr the second)
      % EI + pep --> pyr + EI_P: EI + pep <--> EI_pep <--> EI_P_pyr <--> pyr + EI_P
      % EI_P + HPr --> EI + HPr_P: EI_P + HPr <--> EI_P_HPr <--> EI_HPr_P <--> HPr_P + EI
      % EIIA + HPr_P --> HPr + EIIA_P: EIIA + HPr_P <--> EIIA_HPr_P <--> 
      % EIIBC + EIIA_P --> EIIA + EIIBC_P: EIIA + EIIBC <--> EIIA_P_EIIBC <--> EIIA_EIIBC_P <--> EIIA + EIIBC_P
      % EIIBC_P + glc --> g6p + EIIBC:  EIIBC_P + (glc) <--> EIIBC_P_glc <--> EIIBC_g6p <--> EIIBC + g6p
      elseif [model.rxntype(j)==4] && [max(size(find(model.S(:,j) < 0))) <= 2] && [max(size(find(model.S(:,j) > 0))) == 2]
         reactant_indices = find(model.S(:,j) < 0);
         product_indices = find(model.S(:,j) > 0);

         glc_index = find(strcmp(model.metab,'Glucose')==1);
         g6p_index = find(strcmp(model.metab,'G6P')==1);
         pep_index = find(strcmp(model.metab,'PEP')==1);
         g6p_index = find(strcmp(model.metab,'G6P')==1);
         pyr_index = find(strcmp(model.metab,'Pyruvate')==1);

         EI_index = find(strcmp(model.Enzyme,'EI')==1);
         HPr_index = find(strcmp(model.Enzyme,'HPr')==1);
         EIIA_index = find(strcmp(model.Enzyme,'EIIA')==1);
         EIIBC_index = find(strcmp(model.Enzyme,'EIIBC')==1);
      
 
         if isempty(glc_index) && [max(size(find(model.S(:,j) < 0))) == 2]% By Ali: It's checked if there is just one metabolite then that would be ok, and the second metabolite should be Glucose.
           fprintf('\nError! There is no metabolite named Glucose in the model for the PTS reaction\n');
           return
         elseif isempty(g6p_index)
           fprintf('\nError! There is no metabolite named G6P in the model for the PTS reaction\n');
           return
         elseif isempty(pep_index)
           fprintf('\nError! There is no metabolite named PEP in the model for the PTS reaction\n');
           return
         elseif isempty(g6p_index)
           fprintf('\nError! There is no metabolite named G6P in the model for the PTS reaction\n');
           return
         elseif isempty(pyr_index)
           fprintf('\nError! There is not metabolite named Pyruvate in the model for the PTS reaction\n');
           return
         elseif isempty(EI_index)
           fprintf('\nError! There is no enzyme named EI in the model for the PTS reaction\n');
           return
         elseif isempty(HPr_index)
           fprintf('\nError! There is no enzyme named HPr in the model for the PTS reaction\n');
           return
         elseif isempty(EIIA_index)
           fprintf('\nError! There is no enzyme named EIIA in the model for the PTS reaction\n');
           return
         elseif isempty(EIIBC_index)
           fprintf('\nError! There is no enzyme named EIIBC in the model for the PTS reaction\n');
           return
         end


         if [max(size(find(model.S(:,j) < 0))) == 2] && [glc_index ~+ reactant_indeces(1)] && [glc_index ~= reactant_indices(2)]
           fprintf('\nError! The index of Glucose the fileds metab and S are not consistent \n');
           return
         elseif [pep_index ~= reactant_indices(1)] && [pep_index ~= reactant_indices(2)]
           fprintf('\nError! The index of PEP the fileds metab and S are not consistent \n');
           return
         elseif [g6p_index ~= product_indices(1)] && [g6p_index ~= product_indices(2)]
           fprintf('\nError! The index of G6P the fileds metab and S are not consistent \n');
           return
         elseif [pyr_index ~= product_indices(1)] && [pyr_index ~= product_indices(2)]
           fprintf('\nError! The index of Pyruvate the fileds metab and S are not consistent \n');
           return
         end

         % Create rxn name j_1: EI + pep <--> EI_pep 
         rxn_index = rxn_index + 1;
         expanded_model.rxn(rxn_index) = strcat(model.rxn(j),'_1');

         % Create stoichiomteric matrix: S(pep,J_1)=-1 and S(EI,J_1)=-1
         expanded_model.S(pep_index,rxn_index) = -1;
         expanded_model.S(metab_num_orig+EI_index,rxn_index) = -1;

         expanded_model.rev(rxn_index)=1;
         expanded_model.LB_rxn(rxn_index)=-1000;
         expanded_model.UB_rxn(rxn_index)=1000;
         expanded_model.isElementary(rxn_index)=1;
         expanded_model.unexp_rxn_info(j,2)={rxn_index};  % Start index

         % EI_pep
         metab_index = metab_index + 1;
         expanded_model.metab(metab_index)=strcat(model.Enzyme(EI_index),'_',model.metab(pep_index),'_complex');

         expanded_model.enz_enzComplex(metab_index,EI_index) = 1;
         expanded_model.metab_enzComplex(metab_index,pep_index) = 1;

         % S(EI_pep,j_1)=1
         expanded_model.S(metab_index,rxn_index) = 1;

         % Create rxns J_2 (EI_pep <--> EI_P_pyr) and S(EI_pep,J_2)=-1 
         rxn_index = rxn_index + 1;
         expanded_model.rxn(rxn_index) = strcat(model.rxn(j),'_2');
         expanded_model.S(metab_index,rxn_index) = -1;      

         expanded_model.rev(rxn_index)=1;
         expanded_model.LB_rxn(rxn_index)=-1000;
         expanded_model.UB_rxn(rxn_index)=1000;
         expanded_model.isElementary(rxn_index)=1;
 
         % EI_P_pyr
         metab_index = metab_index + 1;
         expanded_model.metab(metab_index)=strcat(model.Enzyme(EI_index),'_P_',model.metab(pyr_index),'_complex');

         expanded_model.enz_enzComplex(metab_index,EI_index) = 1;
         expanded_model.metab_enzComplex(metab_index,pyr_index) = 1;

         % S(EI_P_pyr,j_2)=1
         expanded_model.S(metab_index,rxn_index) = 1;      

         % Create rxns J_3 (EI_P_pyr <--> pyr + EI_P) and S(EI_P_pyr,J_3)=-1 
         rxn_index = rxn_index + 1;
         expanded_model.rxn(rxn_index) = strcat(model.rxn(j),'_3');
         expanded_model.S(metab_index,rxn_index) = -1;        % S(EI_P_pyr,j_3)=-1

         expanded_model.rev(rxn_index)=1;
         expanded_model.LB_rxn(rxn_index)=-1000;
         expanded_model.UB_rxn(rxn_index)=1000;
         expanded_model.isElementary(rxn_index)=1;

         % EI_P
         metab_index = metab_index + 1;
         expanded_model.metab(metab_index)=strcat(model.Enzyme(EI_index),'_P_complex');

         expanded_model.enz_enzComplex(metab_index,EI_index) = 1;

         % S(EI_P,J_3)=1
         expanded_model.S(metab_index,rxn_index) = 1;

         % S(pyr,J_3)=1
         expanded_model.S(pyr_index,rxn_index) = 1;

         % Create rxn name j_4: EI_P + HPr <--> EI_P_HPr 
         rxn_index = rxn_index + 1;
         expanded_model.rxn(rxn_index) = strcat(model.rxn(j),'_4');

         expanded_model.rev(rxn_index)=1;
         expanded_model.LB_rxn(rxn_index)=-1000;
         expanded_model.UB_rxn(rxn_index)=1000;
         expanded_model.isElementary(rxn_index)=1;

         % S(EI_P,J_4) = -1
         expanded_model.S(metab_index,rxn_index) = -1;

         % S(HPr,J_4) = -1)
         expanded_model.S(metab_num_orig+HPr_index,rxn_index) = -1;

         % EI_P_HPr
         metab_index = metab_index + 1;
         expanded_model.metab(metab_index)=strcat(model.Enzyme(EI_index),'_P_',model.Enzyme(HPr_index),'_complex');

         expanded_model.enz_enzComplex(metab_index,EI_index) = 1;
         expanded_model.enz_enzComplex(metab_index,HPr_index) = 1;

         % S(EI_P_HPr,j_4)=1  
         expanded_model.S(metab_index,rxn_index) = 1;

         % Create rxn name j_5: EI_P_HPr <--> EI_HPr_P
         rxn_index = rxn_index + 1;
         expanded_model.rxn(rxn_index) = strcat(model.rxn(j),'_5');

         expanded_model.rev(rxn_index)=1;
         expanded_model.LB_rxn(rxn_index)=-1000;
         expanded_model.UB_rxn(rxn_index)=1000;
         expanded_model.isElementary(rxn_index)=1;

         % S(EI_P_HPr,j_5)=-1
         expanded_model.S(metab_index,rxn_index) = -1;

         % EI_HPr_P
         metab_index = metab_index + 1;
         expanded_model.metab(metab_index)=strcat(model.Enzyme(EI_index),'_',model.Enzyme(HPr_index),'_P_complex');
         expanded_model.enz_enzComplex(metab_index,EI_index) = 1;
         expanded_model.enz_enzComplex(metab_index,HPr_index) = 1;

         % S(EI_HPr_P,j_5)=1
         expanded_model.S(metab_index,rxn_index) = 1;

         % Create rxn name j_6: EI_HPr_P <--> EI + HPr_P
         rxn_index = rxn_index + 1;
         expanded_model.rxn(rxn_index) = strcat(model.rxn(j),'_6');

         expanded_model.rev(rxn_index)=1;
         expanded_model.LB_rxn(rxn_index)=-1000;
         expanded_model.UB_rxn(rxn_index)=1000;
         expanded_model.isElementary(rxn_index)=1;

         % S(EI_HPr_P,J_6)=-1    and S(EI,j_6)=1     
         expanded_model.S(metab_index,rxn_index) = -1;
         expanded_model.S(metab_num_orig+EI_index,rxn_index) = 1;

         % HPr_P
         metab_index = metab_index + 1;
         expanded_model.metab(metab_index)=strcat(model.Enzyme(HPr_index),'_P_complex');
         expanded_model.enz_enzComplex(metab_index,HPr_index) = 1;

         % S(HPr_P,j_6)=1
         expanded_model.S(metab_index,rxn_index) = 1;

         % Create rxn name j_7: HPr_P + EIIA <--> HPr_P_EIIA
         rxn_index = rxn_index + 1;
         expanded_model.rxn(rxn_index) = strcat(model.rxn(j),'_7');

         expanded_model.rev(rxn_index)=1;
         expanded_model.LB_rxn(rxn_index)=-1000;
         expanded_model.UB_rxn(rxn_index)=1000;
         expanded_model.isElementary(rxn_index)=1;

         % S(HPr_P,J_7)=-1    and S(EIIA,j_7)=-1     
         expanded_model.S(metab_index,rxn_index) = -1;
         expanded_model.S(metab_num_orig+EIIA_index,rxn_index) = -1;

         % HPr_P_EIIA
         metab_index = metab_index + 1;
         expanded_model.metab(metab_index)=strcat(model.Enzyme(HPr_index),'_P_',model.Enzyme(EIIA_index),'_complex');
         expanded_model.enz_enzComplex(metab_index,HPr_index) = 1;
         expanded_model.enz_enzComplex(metab_index,EIIA_index) = 1;

         % S(HPr_P_EIIA,j_7)=1
         expanded_model.S(metab_index,rxn_index) = 1;

         % Create rxn name j_8: HPr_P_EIIA <--> HPr_EIIA_P
         rxn_index = rxn_index + 1;
         expanded_model.rxn(rxn_index) = strcat(model.rxn(j),'_8');

         expanded_model.rev(rxn_index)=1;
         expanded_model.LB_rxn(rxn_index)=-1000;
         expanded_model.UB_rxn(rxn_index)=1000;
         expanded_model.isElementary(rxn_index)=1;

         % S(HPr_P_EIIA,j_8)=-1
         expanded_model.S(metab_index,rxn_index) = -1;

         % HPr_EIIA_P
         metab_index = metab_index + 1;
         expanded_model.metab(metab_index)=strcat(model.Enzyme(HPr_index),'_',model.Enzyme(EIIA_index),'_P_complex');
         expanded_model.enz_enzComplex(metab_index,HPr_index) = 1;
         expanded_model.enz_enzComplex(metab_index,EIIA_index) = 1;

         % S(HPr_EIIA_P,j_8)=1
         expanded_model.S(metab_index,rxn_index) = 1;

         % Create rxn name j_9: HPr_EIIA_P <--> HPr + EIIA_P
         rxn_index = rxn_index + 1;
         expanded_model.rxn(rxn_index) = strcat(model.rxn(j),'_9');

         expanded_model.rev(rxn_index)=1;
         expanded_model.LB_rxn(rxn_index)=-1000;
         expanded_model.UB_rxn(rxn_index)=1000;
         expanded_model.isElementary(rxn_index)=1;

         % S(HPr_EIIA_P,j_9)=-1    S(HPr,J_9)=1
         expanded_model.S(metab_index,rxn_index) = -1;
         expanded_model.S(metab_num_orig+HPr_index,rxn_index) = 1;

         % EIIA_P
         metab_index = metab_index + 1;
         expanded_model.metab(metab_index)=strcat(model.Enzyme(EIIA_index),'_P_complex');
         expanded_model.enz_enzComplex(metab_index,EIIA_index) = 1;

         % S(EIIA_P,J_9) = 1
         expanded_model.S(metab_index,rxn_index) = 1;


         % Create rxn name j_10: EIIA_P + EIIBC <--> EIIA_P_EIIBC
         rxn_index = rxn_index + 1;
         expanded_model.rxn(rxn_index) = strcat(model.rxn(j),'_10');

         expanded_model.rev(rxn_index)=1;
         expanded_model.LB_rxn(rxn_index)=-1000;
         expanded_model.UB_rxn(rxn_index)=1000;
         expanded_model.isElementary(rxn_index)=1;

         % S(EIIA_P,J_10)=-1    and S(EIIBC,j_10)=-1     
         expanded_model.S(metab_index,rxn_index) = -1;
         expanded_model.S(metab_num_orig+EIIBC_index,rxn_index) = -1;

         % EIIA_P_EIIBC
         metab_index = metab_index + 1;
         expanded_model.metab(metab_index)=strcat(model.Enzyme(EIIA_index),'_P_',model.Enzyme(EIIBC_index),'_complex');
         expanded_model.enz_enzComplex(metab_index,EIIA_index) = 1;
         expanded_model.enz_enzComplex(metab_index,EIIBC_index) = 1;

         % S(EIIA_P_EIIBC,j_7)=1
         expanded_model.S(metab_index,rxn_index) = 1;

         % Create rxn name j_11: EIIA_P_EIIBC <--> EIIA_EIIBC_P
         rxn_index = rxn_index + 1;
         expanded_model.rxn(rxn_index) = strcat(model.rxn(j),'_11');

         expanded_model.rev(rxn_index)=1;
         expanded_model.LB_rxn(rxn_index)=-1000;
         expanded_model.UB_rxn(rxn_index)=1000;
         expanded_model.isElementary(rxn_index)=1;

         % S(EIIA_P_EIIBC,j_11)=-1
         expanded_model.S(metab_index,rxn_index) = -1;

         % EIIA_EIIBC_P
         metab_index = metab_index + 1;
         expanded_model.metab(metab_index)=strcat(model.Enzyme(EIIA_index),'_',model.Enzyme(EIIBC_index),'_P_complex');
         expanded_model.enz_enzComplex(metab_index,EIIA_index) = 1;
         expanded_model.enz_enzComplex(metab_index,EIIBC_index) = 1;

         % S(EIIA_EIIBC_P,j_11)=1
         expanded_model.S(metab_index,rxn_index) = 1;

         % Create rxn name j_12: EIIA_EIIBC_P <--> EIIA + EIIBC_P
         rxn_index = rxn_index + 1;
         expanded_model.rxn(rxn_index) = strcat(model.rxn(j),'_12');

         expanded_model.rev(rxn_index)=1;
         expanded_model.LB_rxn(rxn_index)=-1000;
         expanded_model.UB_rxn(rxn_index)=1000;
         expanded_model.isElementary(rxn_index)=1;

         % S(EIIA_EIIBC_P,j_12)=-1    S(EIIA,J_12)=1
         expanded_model.S(metab_index,rxn_index) = -1;
         expanded_model.S(metab_num_orig+EIIA_index,rxn_index) = 1;

         % EIIBC_P
         metab_index = metab_index + 1;
         expanded_model.metab(metab_index)=strcat(model.Enzyme(EIIBC_index),'_P_complex');
         expanded_model.enz_enzComplex(metab_index,EIIBC_index) = 1;

         % S(EIIBC_P,J_9) = 1
         expanded_model.S(metab_index,rxn_index) = 1;

         % Create rxn name j_13: EIIBC_P + glc <--> EIIBC_P_glc
         rxn_index = rxn_index + 1;
         expanded_model.rxn(rxn_index) = strcat(model.rxn(j),'_13');

         expanded_model.rev(rxn_index)=1;
         expanded_model.LB_rxn(rxn_index)=-1000;
         expanded_model.UB_rxn(rxn_index)=1000;
         expanded_model.isElementary(rxn_index)=1;

         % S(EIIBC_P,J_13)=-1    
         expanded_model.S(metab_index,rxn_index) = -1;

         % S(glc,J_13)=-1
         if max(size(find(model.S(:,j) < 0))) == 2
            expanded_model.S(glc_index,rxn_index) = -1;
         end

         % EIIBC_P_glc
         metab_index = metab_index + 1;
         if max(size(find(model.S(:,j) < 0))) == 2
            expanded_model.metab(metab_index)=strcat(model.Enzyme(EIIBC_index),'_P_',model.metab(glc_index),'_complex');
            expanded_model.metab_enzComplex(metab_index,glc_index) = 1;
         elseif max(size(find(model.S(:,j) < 0))) == 1
            expanded_model.metab(metab_index)=strcat(model.Enzyme(EIIBC_index),'_P_','Glc_ext_complex');
         end
         expanded_model.enz_enzComplex(metab_index,EIIA_index) = 1;

         % S(EIIBC_P_glc,j_7)=1
         expanded_model.S(metab_index,rxn_index) = 1;

         % Create rxn name j_14: EIIBC_P_glc <--> EIIBC_G6P
         rxn_index = rxn_index + 1;
         expanded_model.rxn(rxn_index) = strcat(model.rxn(j),'_14');

         expanded_model.rev(rxn_index)=1;
         expanded_model.LB_rxn(rxn_index)=-1000;
         expanded_model.UB_rxn(rxn_index)=1000;
         expanded_model.isElementary(rxn_index)=1;

         % S(EIIBC_P_glc,j_14)=-1
         expanded_model.S(metab_index,rxn_index) = -1;

         % EIIBC_g6p
         metab_index = metab_index + 1;
         expanded_model.metab(metab_index)=strcat(model.Enzyme(EIIBC_index),'_',model.metab(g6p_index),'_complex');
         expanded_model.enz_enzComplex(metab_index,EIIBC_index) = 1;
         expanded_model.metab_enzComplex(metab_index,g6p_index) = 1;

         % S(EIIA_g6p,j_14)=1
         expanded_model.S(metab_index,rxn_index) = 1;

         % Create rxn name j_15: EIIBC_g6p <--> EIIBC + g6p
         rxn_index = rxn_index + 1;
         expanded_model.rxn(rxn_index) = strcat(model.rxn(j),'_15');

         expanded_model.rev(rxn_index)=1;
         expanded_model.LB_rxn(rxn_index)=-1000;
         expanded_model.UB_rxn(rxn_index)=1000;
         expanded_model.isElementary(rxn_index)=1;

         % S(EIIBC_g6p,j_15)=-1   S(EIIBC,J_15)=1     S(g6p,J_15)=1 
         expanded_model.S(metab_index,rxn_index) = -1;
         expanded_model.S(metab_num_orig+EIIBC_index,rxn_index) = 1;
         expanded_model.S(g6p_index,rxn_index) = 1;
         expanded_model.unexp_rxn_info(j,3)={rxn_index}; % end index

      end   % end of if [model.rxntype(j)==1] && ... elseif ...  elseif ... else ...

  % --- Non-enzymatic (Exchange export) rxns -----
  elseif sum(model.GPR(:,j)) == 0   
      reactant_index=find(model.S(:,j) < 0);
      for k=1:max(size(reactant_index)) 
         rxn_index=rxn_index+1;
         expanded_model.rxn(rxn_index) = model.rxn(j); 
         expanded_model.S(reactant_index(k),rxn_index) = model.S(reactant_index(k),j);

         expanded_model.rev(rxn_index)=0;
         expanded_model.LB_rxn(rxn_index)=0;
         expanded_model.UB_rxn(rxn_index)=1000;
         expanded_model.isElementary(rxn_index)=0;
         expanded_model.unexp_rxn_info(j,2)={rxn_index}; % start index
         expanded_model.unexp_rxn_info(j,3)={rxn_index}; % end index
      end    % end of for k=1:max(size(reactant_index))
  end     % end of if sum(model.GPR(:,j)) > 0
end    % end of or j=1:rxn_num_orig

% number of elementary rxns not associated with regulation
noRegulation_rxn_num=rxn_index;

reg_counter=rxn_num_orig;

%-------------- Regulatory reactions ------------ 
for j=1:size(model.enzymeReg,2)        % Loop over enzymes
   for i=1:size(model.enzymeReg,1)     % Loop over metabolites

      %--------------- Competitive inhibition -----------------
      % A + E <--> EA <--> EB <--> E + B
      % I + E <--> EI
      if model.enzymeReg(i,j) == -1 

         reg_counter=reg_counter+1;
         % Original rxn name
         expanded_model.unexp_rxn_info(reg_counter,1)=strcat(model.Enzyme(j),'_comp_inhib_by_',model.metab(i));
         % min deltaG
         expanded_model.unexp_rxn_info(reg_counter,4)={0};
         % max deltaG
         expanded_model.unexp_rxn_info(reg_counter,5)={0};
         % rxntype 
         expanded_model.unexp_rxn_info(reg_counter,6)={-1};
         % steady-state fluxes for the wildtype and perturbed strains
         for np = 1:num_perturb
           expanded_model.unexp_rxn_info(reg_counter,6+np)={0};
         end 

         % Create rxn name enzymeName_comp_inhibition_inhibitorName: E + I <--> E_I
         rxn_index = rxn_index + 1;
         expanded_model.rxn(rxn_index) = strcat(model.Enzyme(j),'_comp_inhibit_by_',model.metab(i));

         expanded_model.rev(rxn_index)=1;
         expanded_model.LB_rxn(rxn_index)=-1000;
         expanded_model.UB_rxn(rxn_index)=1000;
         expanded_model.isElementary(rxn_index)=1;
         expanded_model.unexp_rxn_info(reg_counter,2)={rxn_index};  % Start index
         expanded_model.unexp_rxn_info(reg_counter,3)={rxn_index};  % End index

         % S(E,j) = -1
         expanded_model.S(metab_num_orig + j,rxn_index) = -1;

         % S(I,j) = -1 
         expanded_model.S(i,rxn_index) = -1;

         % EI
         metab_index = metab_index + 1;
         expanded_model.metab(metab_index)=strcat(model.Enzyme(j),'_',model.metab(i),'_complex');

         expanded_model.enz_enzComplex(metab_index,j) = 1;
         expanded_model.metab_enzComplex(metab_index,i) = 1;

         % S(EI,j) = 1
         expanded_model.S(metab_index,rxn_index) = 1;

      %---------------- Uncompetitive inhibition ------------------
      % A + E <--> EA <--> EB <--> E + B
      % I + EA <--> EAI
      elseif model.enzymeReg(i,j) == -2 

        reg_counter + reg_counter + 1;

        % original rxn name 
        expanded_model.unexp_rxn_info(reg_counter,1)=strcat(model.Enzyme(j),'_uncomp_inhibit_by_',expanded_model.metab(i));
        % min deltaG
        expanded_model.unexp_rxn_info(reg_counter,4)={0};
        % max deltaG
        expanded_model.unexp_rxn_info(reg_counter,5)={0};
        % rxntype 
        expanded_model.unexp_rxn_info(reg_counter,6)={-2};
        % steady-state fluxes for the wildtype and perturbed strains
        for np = 1:num_perturb
          expanded_model.unexp_rxn_info(reg_counter,6+np)={0};
        end 

         % Find the index of rxns in which E appears as a reactant
         E_reactant_rxn_index = find(expanded_model.S(j+metab_num_orig,1:noRegulation_rxn_num) < 0);

         for k=1:max(size(E_reactant_rxn_index))
           % Find the index of EA 
           EA_index=find(expanded_model.S(metab_num_orig + enzyme_num+1:end,E_reactant_rxn_index(k)) > 0) + (metab_num_orig + enzyme_num);
          
 
           % Create rxn name enzymeName_uncomp_inhibition_k_inhibitorName: EA + I <--> EAI
           rxn_index = rxn_index + 1;
           expanded_model.rxn(rxn_index) = strcat(model.Enzyme(j),'_uncomp_inhibit_by_',expanded_model.metab(i),'_',num2str(k));

           expanded_model.rev(rxn_index)=1;
           expanded_model.LB_rxn(rxn_index)=-1000;
           expanded_model.UB_rxn(rxn_index)=1000;
           expanded_model.isElementary(rxn_index)=1;
           if k==1
             expanded_model.unexp_rxn_info(reg_counter,2)={rxn_index};  % Start index
           end  

           % S(I,j) = -1
           expanded_model.S(i,rxn_index)=-1;

           % S(EA,j) = -1
           expanded_model.S(EA_index,rxn_index)=-1;

           % EAI: First remove _complex from EA_complex and then add strings together
           metab_index = metab_index + 1;
           expanded_model.metab(metab_index)=strcat(regexprep(expanded_model.metab(EA_index(k)),'_complex',''),'_',model.metab(i),'_complex');

           expanded_model.enz_enzComplex(metab_index,j) = 1;
           expanded_model.metab_enzComplex(metab_index,i) = 1;

           % S(EAI,j)=1
           expanded_model.S(metab_index,rxn_index)=1;

         end  % for k

         % End index
         expanded_model.unexp_rxn_info(reg_counter,3)={rxn_index};  

      %------------------ Mixed inhibition ----------------------
      elseif model.enzymeReg(i,j) == -3 

         reg_counter = reg_counter + 1;

         % Original rxn name
         expanded_model.unexp_rxn_info(reg_counter,1)=strcat(model.Enzyme(j),'_mixed_inhibit_by_',model.metab(i));

         % min deltaG
         expanded_model.unexp_rxn_info(reg_counter,4)={0};
         % max deltaG
         expanded_model.unexp_rxn_info(reg_counter,5)={0};
         % rxntype 
         expanded_model.unexp_rxn_info(reg_counter,6)={-3};
         % steady-state fluxes for the wildtype and perturbed strains
         for np = 1:num_perturb
           expanded_model.unexp_rxn_info(reg_counter,6+np)={0};
         end 

         %-- Competetive part --
         % Create rxn name enzymeName_comp_inhibition_inhibitorName: E + I <--> E_I
         rxn_index = rxn_index + 1;
         expanded_model.rxn(rxn_index) = strcat(model.Enzyme(j),'_comp_inhibit_by_',model.metab(i));

         expanded_model.rev(rxn_index)=1;
         expanded_model.LB_rxn(rxn_index)=-1000;
         expanded_model.UB_rxn(rxn_index)=1000;
         expanded_model.isElementary(rxn_index)=1;
         expanded_model.unexp_rxn_info(reg_counter,2)={rxn_index};  % Start index

         % S(E,j) = -1
         expanded_model.S(metab_num_orig + j,rxn_index) = -1;

         % S(I,j) = -1 
         expanded_model.S(i,rxn_index) = -1;

         % EI
         metab_index = metab_index + 1;
         expanded_model.metab(metab_index)=strcat(model.Enzyme(j),'_',model.metab(i),'_complex');

         expanded_model.enz_enzComplex(metab_index,j) = 1;
         expanded_model.metab_enzComplex(metab_index,i) = 1;

         % S(EI,j) = 1
         expanded_model.S(metab_index,rxn_index) = 1;

         %-- Uncompetetive part --
         % Find the index of rxns in which E appears as a reactant

         E_reactant_rxn_index = find(expanded_model.S(j+metab_num_orig,1:noRegulation_rxn_num) < 0);

         for k=1:max(size(E_reactant_rxn_index))
           % Find the index of EA 
           EA_index=find(expanded_model.S(metab_num_orig + enzyme_num+1:end,E_reactant_rxn_index(k)) > 0) + (metab_num_orig + enzyme_num);
      
           % Create rxn name enzymeName_uncomp_inhibition_k_inhibitorName: EA + I <--> EAI
           rxn_index = rxn_index + 1;
           expanded_model.rxn(rxn_index) = strcat(model.Enzyme(j),'_uncomp_inhibit_by_',expanded_model.metab(i),'_',num2str(k));

           expanded_model.rev(rxn_index)=1;
           expanded_model.LB_rxn(rxn_index)=-1000;
           expanded_model.UB_rxn(rxn_index)=1000;
           expanded_model.isElementary(rxn_index)=1;

           % S(I,j) = -1
           expanded_model.S(i,rxn_index)=-1;

           % S(EA,j) = -1
           expanded_model.S(EA_index,rxn_index)=-1;

           % EAI: First remove _complex from EA_complex and then add strings together
           metab_index = metab_index + 1;
           expanded_model.metab(metab_index)=strcat(regexprep(expanded_model.metab(EA_index),'_complex',''),'_',model.metab(i),'_complex');

           expanded_model.enz_enzComplex(metab_index,j) = 1;
           expanded_model.metab_enzComplex(metab_index,i) = 1;

           % S(EAI,j)=1
           expanded_model.S(metab_index,rxn_index)=1;

         end  % for k

         % End index
         expanded_model.unexp_rxn_info(reg_counter,3)={rxn_index};  

      % Activation
      elseif model.enzymeReg(i,j) == -4 

        reg_counter = reg_counter + 1;

        % original rxn name 
        expanded_model.unexp_rxn_info(reg_counter,1)=strcat(model.Enzyme(j),'_activ_by_',expanded_model.metab(i));
        % min deltaG
        expanded_model.unexp_rxn_info(reg_counter,4)={0};
        % max deltaG
        expanded_model.unexp_rxn_info(reg_counter,5)={0};
        % rxntype 
        expanded_model.unexp_rxn_info(reg_counter,6)={-4};
        % steady-state fluxes for the wildtype and perturbed strains
        for np = 1:num_perturb
          expanded_model.unexp_rxn_info(reg_counter,6+np)={0};
        end
      end    % end of if model.enzymeReg(i,j) == -1 elseif ... elseif ... else ...
   end    % end of for i=1:size(model.enzymeReg,1) 
end   % end of for j=1:size(model.enzymeReg,2) 

expanded_model.metabIndexInfo(3,:) = {'Enzyme complex' metab_num_orig+enzyme_num+1 max(size(expanded_model.metab))};


expanded_model.rxn=expanded_model.rxn';
expanded_model.rev=expanded_model.rev';
expanded_model.LB_rxn=expanded_model.LB_rxn';
expanded_model.UB_rxn=expanded_model.UB_rxn';
expanded_model.isElementary=expanded_model.isElementary';
expanded_model.Enzyme = model.Enzyme;

% Rxn names for reversible rxns decomposed into forward and backward
rev_num = 0;     % number of reversible rxns
irrev_num=0;     % number of irreversible rxns

for j=1:size(expanded_model.S,2)
     if expanded_model.rev == 1
        rev_num=rev_num+1;
     else
        irrev_num=irrev_num+1;
     end
end


% This field shows what the corresponding index of a reaction in expanded_model.rxn is
% in the expanded_model.rxn_f_b 
expanded_model.rxn_f_bInd_rxnInd = zeros(2*rev_num + irrev_num,1);

expanded_model.rxn_f_b = cell(2*rev_num + irrev_num,1);
counter=0;
for j=1:size(expanded_model.S,2)
  if expanded_model.rev(j) ==1      % if reversible
     counter=counter+1;
     expanded_model.rxn_f_b(counter) = strcat(expanded_model.rxn(j),'_f');
     expanded_model.rxn_f_bInd_rxnInd(counter) = j;
     counter=counter+1;
     expanded_model.rxn_f_b(counter) = strcat(expanded_model.rxn(j),'_b');
     expanded_model.rxn_f_bInd_rxnInd(counter) = j;
  else                             % if irreversible
     counter=counter+1;
     expanded_model.rxn_f_b(counter) = expanded_model.rxn(j);
     expanded_model.rxn_f_bInd_rxnInd(counter) = j;
  end
end


% S and rxn names for reversible rxns decomposed into forward and backward
expanded_model.S_f_b = zeros(size(expanded_model.S,1),2*rev_num + irrev_num);

for i=1:size(expanded_model.S,1)
   counter = 0;
   for j=1:size(expanded_model.S,2)
      if expanded_model.rev(j) ==1      % if reversible
         counter=counter+1;
         expanded_model.S_f_b(i,counter)=expanded_model.S(i,j);;
         counter=counter+1;
         expanded_model.S_f_b(i,counter)=-expanded_model.S(i,j);;
      else                             % if irreversible
         counter=counter+1;
         expanded_model.S_f_b(i,counter)=expanded_model.S(i,j);;
      end
   end
end

% Index of cofactors (the first cofactor only)
expanded_model.cof_metab_index=model.cof_metab_index;

% Index of the rxns converting cofactors (e.g., ATP <--> ADP)
expanded_model.cof_rxn_index = model.cof_rxn_index;

% Perturbations
expanded_model.Perturbations = model.Perturbations;

% Model name
expanded_model.model_name = model.model_name;

% The whole unexpanded model
expanded_model.unexpanded_model = model;

expanded_model=orderfields(expanded_model,[1 5 13 6 16 7 8 9 10 15 3 4 2 17 18 12 14 19 11 20 21]);
