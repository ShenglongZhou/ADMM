function [Arank] = Approx_rank(A)
     p    = size(A,1);
    [~,D] = svd(A); 
    for k = 2:p
        if D(k,k)/D(1,1)<0.01; break; end
    end
    if k == p;  Arank=p; else; Arank=k-1; end
end