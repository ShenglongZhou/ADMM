function [Sigma, time] = ADMM(Sigman,opts)
% A solver for sparse and low rank covariance matrix recovery model:
%
% min 0.5||Sigma-Sigman||_F^2 + lambda*||Sigma||_1+tau*||Sigma||_*
%
% Inputs:
%    Sigman --- an p x p matrix; (Required)
%    opts   --- a structure with several fields: (Optional)
%        opts.lambda -- penalty for sparsity, default: 0.5
%        opts.tau    -- penalty for low rank, default: 0.5                               
%        opts.tol    -- tolerance for halting condition, default: 5e-4 
%        opts.itmax  -- maximum number of iterations, default: 1000 
%
% Outputs:
%    Sigma  --- an p x p matrix, recovered covariance matrix
%    time   --- CPU time
%
% Last updated by Shenglong Zhou, 05/07/2020. 
% This Matlab solver was created based on the algorithm proposed by
% S. Zhou, N. Xiu, Z. Luo and L. Kong, (2015),
% Sparse and Low-Rank Covariance Matrix Estimation, 
% Journal of the Operations Research Society of China, 3(2): 231-250.
% Send your comments and suggestions to  <<<longnan_zsl@163.com >>>                                      %%%%%% 
% Warning: Accuracy may not be guaranteed !!!!!  
% -------------------------------------------------------------------------

% Initialization
if nargin ==1; opts = []; end
[lambda,tau,tol,itmax,mu,mu1,rate,nT] = Parameters();
if isfield(opts,'lambda'); lambda = opts.lambda;  end
if isfield(opts,'tau');    tau    = opts.tau;     end 
if isfield(opts,'tol');    tol    = opts.tol;     end 
if isfield(opts,'itmax');  itmax  = opts.itmax;   end 

p        = size(Sigman,2);   
Lamda    = zeros(p); 
Sigma    = zeros(p);  
Gamma    = speye(p); 
ErrGamma = Inf; 
mutau    = mu*tau;
Fnorm    = @(x)(norm(x,'fro'));
to       = tic; 
fprintf(' Iter     CPUTime     ErrSigma     ErrGamma  \n');
fprintf('--------------------------------------------\n');             
for iter=1:itmax

    % Update Gamma
    Gamma1 = Gamma;
    apsp   = ceil(10+rate*nT); 
    if iter> 1 &&  apsp < ceil(0.2*p)   
    [U,Q]  = eigs(Sigma+mu*Lamda-mutau*speye(p),apsp,'largestreal');
    else
    [U,Q]  = eig(Sigma+mu*Lamda-mutau*speye(p)); 
    end
    Q      = real(Q);
    if size(Q,2)==1
    dQ     = Q;
    else
    dQ     = diag(Q);
    end    
    
    U      = real(U); 
    T      = find(dQ>1e-4);
    nT     = nnz(T);
    Gamma  = (U(:,T).*repmat(dQ(T)',[p,1]))*U(:,T)'; 

    % Update Sigma
    P      = Sigman + Gamma/mu - Lamda; 
    Sigma1 = Sigma;   
    Sigma  = sign(P).* max(abs(P)*mu1-lambda*mu1,0);

    % Compute the Stop Criteria 
    time      = toc(to);
    ErrGamma0 = ErrGamma;
    ErrSigma  = Fnorm(Sigma-Sigma1)/(1+Fnorm(Sigma1));    
    ErrGamma  = Fnorm(Gamma1-Gamma)/(1+Fnorm(Gamma1));      
    fprintf(' %3d      %5.2f       %6.2e     %6.2e\n',...
                         iter,time,ErrSigma,ErrGamma);
    if  ErrSigma<tol && ErrGamma<tol; break;  end    
    
    if  ErrGamma0<ErrGamma  
        rate = rate+0.1; 
    else
        rate = 1;
    end
    
    % Update Lamda
    Lamda = Lamda - (Gamma-Sigma)/mu;
    clear U Q P Sigma1 Gamma1 
end
fprintf('--------------------------------------------\n');     
end

% set parameters ----------------------------------------------------------
function [lamda,tau,tol,itmax,mu,mu1,rate,nT]=Parameters() 
    lamda = 0.5;
    tau   = 0.5;
    tol   = 5e-4;
    itmax = 1000;
    mu    = 1;
    mu1   = mu/(mu+1);
    rate  = 1;
    nT    = 0;
end
