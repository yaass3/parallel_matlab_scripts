function out=brute
    global best_y iteration
    clc 
    delete intermediateGA.mat WholeGA.mat
    best_y=10000;
    iteration=0;

    if matlabpool('size') == 0 % checking to see if my pool is already open
    matlabpool open 2
    end
% tic;
    parfor i=1:2   %28
        
        for h=i+1:2
            for w=h+1:3
        k=1;
        for j=[1e-7,2]
            kk=1;
            for z=1e-7%[1e-7,2]
                kkk=1;
                for r=1e-7%[1e-7,2]
                
            Enzyme_IDs=[i,h,w];
            perturbation_values=[j,z,r];
            
%             iteration=iteration+1;
            A = EM_mainModule_hypothetical_base(Enzyme_IDs,perturbation_values)
%          toc
%              y(i,h,w,k,kk,kkk)= A
    min_obj=A; %-y(i,h,w,k,kk,kkk);
%     iteration
    if min_obj< best_y
       disp 'yesss!!!!!!!!!!!   update yessssssssssss'
       
       best_y=min_obj;
       best_Enzyme_ID=Enzyme_IDs;
       best_perturbation_values=perturbation_values;
       save('intermediateGA.mat','best_y','best_Enzyme_ID','best_perturbation_values'); 
      
    end
    kkk=kkk+1;
                end
    
    kk=kk+1;
            end
            
            k=k+1;
        end
            
            disp '*************************************'
        end
    end
%     save('WholeGA.mat','y','iteration');
    end
%toc;
matlabpool close