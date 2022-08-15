% Inverse Power Method

function [m,u,k] = Inverse_Power(A,eps,Nmax,p)
if nargin < 5
    p=0;
end
if nargin < 4
    Nmax=500;
end
if nargin < 3
    eps=1e-5;
end
n=length(A);
A=A-p*eye(n);
u=ones(n,1);
k=0;
m1=0;
Ainv=inv(A);
while k<Nmax
    v=Ainv*u;
    [~,i]=max(abs(v));
    m=v(i);
    u1=u;
    u=v/m;
    if norm(u-u1)<eps
        break;
    end
    m1=m;
    k=k+1;
end
m=1/m+p;
end