% Gauss Elimination Method

function [x] = Gauss_Elimination(a,b)
n=length(b);
x=zeros(n,1);
if rank(a)==n
        for k=1:n-1
            [~,q]=max(abs(a(k:n,k)));
            q=q+k-1;
            t=a(q,:);a(q,:)=a(k,:);a(k,:)=t;
            t=b(q);b(q)=b(k);b(k)=t;
            for j=(k+1):n
                s=a(j,k)/a(k,k);
                a(j,k:n)=a(j,k:n)-a(k,k:n)*s;
                b(j)=-b(k)*s+b(j);
            end
        end
        x(n)=b(n)/a(n,n);
        for i=(n-1):-1:1
            x(i)=(b(i)-a(i,(i+1):n)*x((i+1):n))/a(i,i);
        end
    else
        disp('rank(a)~=n')
end
end