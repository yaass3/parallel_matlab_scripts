function [jacobFun_sym,f_sym,v_sym,y_vec_sym] = mass_balance_symbolic(expanded_model)
%
% This function computes the RHS of mass balance equations i.e.
% dy/dt = sum(j,S(i,j)*v(j))
%
% INPUTS:
% -------
%   expanded_model:  Expanded model
%
% OUTPUTS:
% ---------
%         jacobFun:  Jacobian matrix
%                f:  Symboic mass balance equations
%                v:  Symbolic elementary reaction rates
%            y_vec:  Symboic y variables       
%
% Ali R. Zomorrodi April 2012
%

S_f_b = expanded_model.S_f_b;

y_vec_sym = sym('y_%d__',[1 length(expanded_model.metab)]);
k_vec_sym = sym('k_%d__',[1 length(expanded_model.rxn_f_b)]);

digits(8);

% --------------------------------------------------
%fprintf('\nComputing rates ....\n');

% First compute the rate of reactions 
for j = 1:length(expanded_model.rxn_f_b)

    v_sym(j) = k_vec_sym(j);

    for i=1:length(expanded_model.metab) 
        if S_f_b(i,j) < 0
            v_sym(j)=v_sym(j)*y_vec_sym(i); 
        end
    end
end

% --------------------------------------------------
%fprintf('\nComputing f_sym ....\n');

for i = 1:length(expanded_model.metab)
  counter = 0;
  for j=1:length(expanded_model.rxn_f_b)
    if S_f_b(i,j) ~= 0
       counter = counter + 1;
       if counter==1
          f_sym(i) = sym(S_f_b(i,j),'d')*v_sym(j);
       else
          f_sym(i) = f_sym(i) + sym(S_f_b(i,j),'d')*v_sym(j);
       end
    end
  end
end


% --------------------------------------------------
%fprintf('\nComputing jacobian ....\n');
jacobFun_sym=jacobian(f_sym,y_vec_sym);

% --------------------------------------------------
%fprintf('\nWriting the results into the file ....\n');
f1=fopen(strcat('jacobFun_',expanded_model.model_name,'.m'),'w');

fprintf(f1,strcat('function df_dy = jacobFun_',expanded_model.model_name,'(t,y,k)\n\n'));
fprintf(f1,'%% This function is created by mass_balance_symbolic.m\n\n');
fprintf(f1,'df_dy = zeros(length(y));\n\n');

for i=1:size(jacobFun_sym,1)
  for j=1:size(jacobFun_sym,2)
     if jacobFun_sym(i,j) ~= 0
        derivative = regexprep(char(jacobFun_sym(i,j)),'__',')');
        derivative = regexprep(derivative,'_','(');
        fprintf(f1,'df_dy(%i,%i) = %s;\n',i,j,derivative);
     end
  end
end
