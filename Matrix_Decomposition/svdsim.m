function [u,s,v] = svdsim(a,tol)
% SVDSIM  simple SVD program
%
% A simple program that demonstrates how to use the
% QR decomposition to perform the SVD of a matrix.
% A may be rectangular and complex.

if ~exist('tol','var')
   tol=eps*1024;
end
   
%reserve space in advance
sizea=size(a);
loopmax=100*max(sizea);
loopcount=0;
% or use Bidiag(A) to initialize U, S, and V
u=eye(sizea(1));
s=a';
v=eye(sizea(2));
Err=realmax;
while Err>tol && loopcount<loopmax
%   log10([Err tol loopcount loopmax]); pause
    [q,s]=qr(s'); u=u*q;
    [q,s]=qr(s'); v=v*q;
% exit when we get "close"
    e=triu(s,1);
    E=norm(e(:));
    F=norm(diag(s));
    if F==0, F=1;end
    Err=E/F;
    loopcount=loopcount+1;
end
% [Err/tol loopcount/loopmax]
% fix the signs in S
ss=diag(s);
s=zeros(sizea);
for n=1:length(ss)
    ssn=ss(n);
    s(n,n)=abs(ssn);
    if ssn<0
       u(:,n)=-u(:,n);
    end
end
if nargout<=1
   u=diag(s);
end
return