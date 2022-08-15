% LU Factorization

function [x]=LU_Factor(A,d)
A
[n,n]=size(A);
L=zeros(n,n);
U=zeros(n,n);
for i=1:n
    L(i,i)=1;
end
for k=1:n
    for j=k:n
        U(k,j)=A(k,j)-sum(L(k,1:k-1).*U(1:k-1,j)');
    end
    for i=k+1:n
        L(i,k)=(A(i,k)-sum(L(i,1:k-1).*U(1:k-1,k)'))/U(k,k);
    end
end

y(1)=d(1);
for i=2:n    
    for j=1:i-1
        d(i)=d(i)-L(i,j)*y(j);
    end
    y(i)=d(i);
end
 
x(n)=y(n)/U(n,n);
for i=(n-1):-1:1
    for j=n:-1:i+1
        y(i)=y(i)-U(i,j)*x(j);
    end
    x(i)=y(i)/U(i,i);
end