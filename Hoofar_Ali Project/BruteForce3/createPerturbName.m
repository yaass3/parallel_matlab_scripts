function perturb_names = createPerturbName(model)

% This function takes a model as input and creates a cell array containing the
% name of pertrurbations
%
% INPUTS:
% -------
%     model = A matlab structure containing the metabolic model
%
% OUTPUTs:
% --------
%     perturb_name =  A cell array containing the name of perturbations
%
% Ali R. Zomorrodi- Costas Maranas Lab, April 2012
%

perturb_names = cell(1,size(model.Perturbations,2));
perturb_names(1)={'wild'};

for j=2:size(model.Perturbations,2)
    perturebed_genes = find(model.Perturbations(:,j) ~= 1);
    if max(size(perturebed_genes)) == 1
       if model.Perturbations(perturebed_genes,j) == 0
         perturb_names(j) = strcat(model.Enzyme(perturebed_genes),'_ko'); 
       elseif model.Perturbations(perturebed_genes,j) < 1
         perturb_names(j) = strcat(model.Enzyme(perturebed_genes),'_dr'); 
       elseif model.Perturbations(perturebed_genes,j) > 1
         perturb_names(j) = strcat(model.Enzyme(perturebed_genes),'_oe'); 
       end
    elseif max(size(perturebed_genes)) > 1
       for k = 1:max(size(perturebed_genes))
          if k==1 
            if model.Perturbations(perturebed_genes(k),j) == 0
               perturb_name = strcat(model.Enzyme(perturebed_genes(k)),'_ko');
            elseif model.Perturbations(perturebed_genes(k),j) < 1
               perturb_name = strcat(model.Enzyme(perturebed_genes(k)),'_dr');
            elseif model.Perturbations(perturebed_genes(k),j) > 1
               perturb_name = strcat(model.Enzyme(perturebed_genes(k)),'_oe');
            end
             
          else  
            if model.Perturbations(perturebed_genes(k),j) == 0
               perturb_name = strcat(perturb_name,'_',model.Enzyme(perturebed_genes(k)),'_ko');
            elseif model.Perturbations(perturebed_genes(k),j) < 1
               perturb_name = strcat(perturb_name,'_',model.Enzyme(perturebed_genes(k)),'_dr');
            elseif model.Perturbations(perturebed_genes(k),j) > 1
               perturb_name = strcat(perturb_name,'_',model.Enzyme(perturebed_genes(k)),'_oe');
            end

          end
       end
       perturb_names(j)=perturb_name;
   end
end


