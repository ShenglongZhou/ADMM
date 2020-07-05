function SubGraph(A,a)

B   = ones(size(A));
absA= abs(A);
T   = find(absA>1e-4);
B(T)= (max(absA(:))-absA(T))/1.75;   

imshow(B)
set(gca,'FontName','Times','FontSize',10)
switch a
    case 1;    title('Population Covariance')
    case 2;    title('Sample Covariance')
    otherwise; title('Recovered Covariance')
end

end