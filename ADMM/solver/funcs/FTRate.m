function [FPR,TPR]=FTRate(A0,A)
 
eps   = 1e-4;
Index = find(abs(A0)>eps);
FPR   = length(find(abs(A(Index))<=eps))/length(find(abs(A)<=eps));
TPR   = length(find(abs(A(Index))>eps))/length(find(abs(A)>eps));

end 