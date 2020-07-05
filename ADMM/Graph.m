clc;clear all;close all;
addpath(genpath(pwd)); warning off

ExType = 1;
p      = 500; 
n      = 50;  
blocks = 5; 

if  ExType==1
    opts.lambda=0.25; opts.tau=0.50;
else
	opts.lambda=0.50; opts.tau=0.75;
end

% Design Sigma0 and Sigman (Sigma0 with Noise) 
[Sigman,Sigma0]= Examples(ExType,n,p,blocks); 

% Call ADMM solver to solve
[Sigma,time]= ADMM(Sigman,opts);

% Result Output
r0        = Approx_rank(Sigma0); 
r         = Approx_rank(Sigma);
Tol       = 1e-4;
sp0       = sum(sum(abs(Sigma0)>=Tol) )/p^2;
sp        = sum(sum(abs(Sigma)>=Tol) )/p^2;
[FPR,TPR] = FTRate(Sigma0,Sigma);
fprintf('CPU Time:  %5.3f(sec); \n',time)
fprintf('AppRankSigma0:  %5d;  AppRankSigma:  %5d\n',r0,r)
fprintf('SparsitySigma0: %5.3f;  SparsitySigma: %5.3f\n',sp0,sp)
fprintf('FalsePosiRate:  %5.3f;  TruePosiRate:  %5.3f\n',FPR,TPR)

% Graph Output
fprintf('--------------------------------------------\n');             
fprintf('Waiting for Graph generation...  \n')
subplot(1,3,1), SubGraph(Sigma0,1);
subplot(1,3,2), SubGraph(Sigman,2);
subplot(1,3,3), SubGraph(Sigma,3);
colormap(summer)


