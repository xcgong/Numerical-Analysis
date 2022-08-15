% Power Method

function [m,x,k] = Power(A,x,eps,Nmax,p)
if nargin < 5
    p=0;
end
if nargin < 4
    Nmax=500;
end
if nargin < 3
    eps=1e-5;
end
n=size(A,2);
if nargin < 2
    x=ones(n,1);
end
A=A-p*eye(n);
k=0;
while k<Nmax
    [~,i]=max(abs(x));
    m=x(i);
    x=x/m;
    x=A*x;
    if k>0
        if abs(m-m1)<eps
            break;
        end
    end
    m1=m;
    k=k+1;
end
m=m+p;
end