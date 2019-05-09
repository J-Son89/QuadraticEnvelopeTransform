function m = Prox2(k,gamma,rho,y)
y0=y;
w=sign(y);
[y1,I]=sort(abs(y),'descend');

Z=zeros(length(y));
for h=1:length(y)
    Z(h,I(h))=1;
end

if (1/(gamma/rho))*y1(k+1)<y1(k)
    m=w.*(Z'*sort(WeightSort(k,rho,gamma,y1),'descend')')';
    m=(rho*y0-gamma*m)/(rho-gamma);

else %We compute the indices called l* and j* in the paper; we call them just l and j.
    y=WeightSort(k,rho,gamma,y1);
    j=1;
    while y(j)>y(k+1)
        j=j+1;
    end

    l=k;
    while  y(l+1)>=y(k)
        l=l+1;
        if l==length(y)
            break
        end
    end

    z=sort(y(j:l),'descend'); %We compute the vector z according to point 3.
    for i=1:length(z)-1
        s=(z(i)+z(i+1))/2;
        x=Candidate(k,s,y);
        j1=1; %we compute the indices called l and j in the paper. We call them l1 and j1.
        while x(j1)~=s
            j1=j1+1;
        end
            l1=length(y);
        while x(l1)~=s
            l1=l1-1;
        end

        sI=(rho*sum(y1(j1:l1)))/((k+1-j1)*rho + (l1-k)*gamma);
        if (sI>=z(i+1) && z(i)>=sI)
                m=w.*(Z'*sort(Candidate(k,sI,y),'descend')')';

            break
        end
    end
    m=(rho*y0-gamma*m)/(rho-gamma);
end