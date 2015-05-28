function E_all = GPHMM_entropy(p)
% p: matrix with the sum of each cloumn is one or 1xN cell array
E_all = [];
    
if ~iscell(p)
    for j = 1:size(p,2)
        E = 0;
        for i=1:size(p,1)
            if (p(i,j)>0)
                E = E - p(i,j)*log2(p(i,j));
            end
        end
        E_all = [E_all E];
    end
end

