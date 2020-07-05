function [Sigman,Sigma0] = Examples(Type,n,p,K)

Sigma0 = zeros(p,p); 
Sigman = zeros(p,p);

if  Type==1
    b1    = zeros(1,K); 
    b1(1) = 1;
    block = randi([floor(p/K/2),floor(3*p/K/2)],1,K-1); 
    if  sum(block)>p  
        block=floor(block-(sum(block)-p)/(K-1)-1); 
    end
    b   = [block,p-sum(block)]';
    for i= 2:K  
        b1(i)= sum(b(1:(i-1)))+1; 
    end
    for i = 1:K
        v =  unifrnd(-1,1,b(i),1);
        Sigma0(b1(i):(b1(i)+b(i)-1),b1(i):(b1(i)+b(i)-1))=v*v';
    end    
elseif Type==2  
    for i = 1:p
    	Sigma0(i,1:p) = max(1-abs(i-(1:p))/10,0);        % %Banding:(1-|i-j|,0)+
    end    
end

Sigma00 = real(Sigma0^(1/2));
X       = Sigma00*randn(p,n);
%X       = mvnrnd(zeros(p,1),Sigma0,n)';
MX      = sum(X,2)/n;
for i=1:n; Sigman = Sigman+(X(:,i)-MX)*(X(:,i)-MX)'; end
Sigman  = Sigman/(n-1);

end


