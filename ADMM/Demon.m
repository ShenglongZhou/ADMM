clc;clear all;close all;
addpath(genpath(pwd)); warning off

ExType = 1;
p      = 500; 
n      = 50;
blocks = 5;   
noSamp = 10;  
result = []; 
 
if  ExType==1
    opts.lambda=0.25; opts.tau=0.50;
else
    opts.lambda=0.50; opts.tau=0.75;
end

for S=1:noSamp
   % Design Sigma0 and Sigman (Sigma0 with Noise) 
    [Sigman,Sigma0]= Examples(ExType,n,p,blocks);
  
    % Call ADMM solver to solve
    [Sigma,time]   = ADMM(Sigman,opts);
    
    r0       = Approx_rank(Sigma0);  
    r        = Approx_rank(Sigma);  
    sp0      = sum(sum(abs(Sigma0)>=1e-4))/p^2;
    sp       = sum(sum(abs(Sigma)>=1e-4))/p^2;
    [FPR,TPR]= FTRate(Sigma0,Sigma);
    
    result   = [result; r0  r sp0 sp FPR TPR  time];
end

% Average Result Output
aver = mean(result,1);
fprintf('\n Rank_S0:  %6.2f     Rank_S: %6.2f',aver(1),aver(2))
fprintf('\n Spar_S0:  %6.4f     Spar_S: %6.4f',aver(3),aver(4))
fprintf('\n FPR:      %6.4f     TPR:    %6.4f',aver(5),aver(6))
fprintf('\n CPU Time: %6.4fsec\n\n',aver(7))
