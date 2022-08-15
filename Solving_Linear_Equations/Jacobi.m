% Jacobi Iteration

function [x, k, index] = Jacobi(A, b, ep, it_max)
if nargin < 5
    it_max=100; 
end
if nargin < 4
    ep=1e-5; 
end
if nargin < 3
    error('Inadequacy of parameters.');
end
n=length(A); 
k=0;
x=zeros(n,1); y=zeros(n,1); 
index=1;
while 1   
    for i=1:n              
        y(i)=b(i);       
        for j=1:n           
            if j~=i               
                y(i)=y(i)-A(i,j)*x(j);           
            end
        end
        if abs(A(i,i))<1e-10 || k==it_max           
            index=0; 
            return;       
        end
        y(i)=y(i)/A(i,i);  
    end
    if norm(y-x,inf)<ep       
        break;   
    end
    k=k+1;
    x=y; 
end
fprintf('%d',k);
end