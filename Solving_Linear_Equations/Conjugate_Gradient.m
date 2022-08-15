% Conjugate Gradient Method

function [x,index,k] = Conjugate_Gradient(A,b,ep,it_max)
k=0;
index=1;
n=size(A,1);
x=zeros(n,1);
r=b-A*x;
p=b-A*x;
while 1
    alpha=(norm(r,2))^2/((A*p)'*p);
    x=x+alpha*p;
    rn=r-alpha*A*p;
    if norm(r,2)<ep
        return
        beta=(norm(rn,2)/norm(r,2))^2;
        r=rn;
        p=r+beta*p;
        return
    end
    beta=(norm(rn,2)/norm(r,2))^2;
    r=rn;
    p=r+beta*p;
    k=k+1;
    if k>it_max
        index=0;
        return
    end
    fprintf('%d',k);
end
end