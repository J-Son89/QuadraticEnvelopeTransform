function F = formSignal(x,C,Zeta,Rank)
    for i = 1:length(x)
        sum=0;
        for j = 1:Rank
            sum = sum +  C(j) * exp( Zeta(j)*x(i) );
        F(i)=sum;
        end
    end
end