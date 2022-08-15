% Jacobi's method for calculation of the eigenvalues and eigenvectors of a symmetric matrix

function [A,Q] =Jacobi_Eigenvalue(A)
n=size(A,1);
ep=1e-5;
tol=ep*n^2;
Q=eye(n);
err=sum(sum(A.^2))-trace(A.^2);
while(err>tol)
    [p,q]=find(abs(A-diag(diag(A)))==max(max(abs(A-diag(diag(A))))),1);
    if abs(A(q,q)-A(p,p))<ep
        theta=sign(A(p,q))*pi/4;
    else
        theta=atan(2*A(p,q)/(A(p,p)-A(q,q)))/2;
    end
    if(abs(theta)>pi/4)
        theta=theta-sign(theta)*pi/2;
    end
    R=[cos(theta) -sin(theta);sin(theta) cos(theta)];
    Q(:,[p,q])=Q(:,[p,q])*R;
    A([p,q],:)=R'*A([p,q],:);
    A(:,[p,q])=A(:,[p,q])*R;
    err=sum(sum(A.^2))-trace(A.^2);
end
end