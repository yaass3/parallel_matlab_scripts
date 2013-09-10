function dy = mass_balance_ode(t,y,kinetic_param,S_f_b)
%
% This function computes the RHS of mass balance equations i.e.
% dy/dt = sum(j,S(i,j)*v(j))
%
% INPUTS:
% -------
%             y: Initial concentrations of metabolites
% kinetic_param: Kinetic values for rxns 
%         S_f_b: Stoichiometric matrix
%
% OUTPUTS:
% ---------
%            dy: Value of dy/dt
%
% Ali R. Zomorrodi April 2012
%


dy = zeros(length(y),1);

% Reaction rates
v = zeros(length(kinetic_param),1);


% First compute the rate of reactions 
for j = 1:length(kinetic_param)

    v(j) = kinetic_param(j);

    for i=1:length(y) 
        if S_f_b(i,j) < 0
            v(j)=v(j)*y(i); 
        end
    end
end


% Loop over each metabolite
for i = 1:length(y)
    dy(i) = S_f_b(i,:)*v;
end

