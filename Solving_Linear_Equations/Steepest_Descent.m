% Steepest Descent Method

function [x,index,k] = Steepest_Descent(A,b,ep,it_max)
k=0;
index=1;
n=size(A,1);
x=zeros(n,1);
r=b-A*x;
p=b-A*x;
while 1
    a=(r'*p)/((A*p)'*p);
    x=x+a*p;
    r=b-A*x;
    p=b-A*x;
    if norm(r,2)<ep
        return
    end
    k=k+1;
    if k>it_max
        index=0;
        return
    end
    fprintf('%d',k);
end
end