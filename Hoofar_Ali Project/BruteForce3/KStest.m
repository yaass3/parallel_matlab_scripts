f=0;
for i=1:100
    for j=1:100
        for k=1:2
            for l=1:28
                if hh(i,j,k,l)<0.01
                    f=f+1;
                end
            end
        end
    end
end
f