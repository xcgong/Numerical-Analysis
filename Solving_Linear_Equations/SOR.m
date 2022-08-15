% SOR Iteration

function [x, k, index] = SOR(A, b, w, ep, it_max)
if nargin < 5
    it_max=100;
end
if nargin < 4
    ep=1e-5;
end
if nargin < 3
    error('Inadequacy of parameters.');
end
[m,n]=size(A);
lenb=length(b);
if m~=n||m~=lenb
    error('The dimension of the equation is incorrect.');
end
n=length(A);
k=0;
x=zeros(n,1);
index=1;
while 1
    x_old=x;
    for i=1:n
        z=b(i);
        for j=1:n
            if j~=i
            z=z-A(i,j)*x_old(j);
        end
    end
    if abs(A(i,i))<1e-10||k==it_max
        index=0;
        return;
    end
    z=z/A(i,i);
    x(i)=z;
    end
    x=w*x+(1-w)*x_old;
    if norm(x_old-x,inf)<ep
        break;
    end
    k=k+1;
end
fprintf('%d',k);
end