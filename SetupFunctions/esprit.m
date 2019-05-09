function Zeta = esprit(Hf,Rank,N)
    Zeta = 1:Rank;
        [U,S,V] = svd(Hf);
        Uplus = U( (2:end) , 1:Rank);
        Uminus =U( (1:end-1) , 1:Rank);
        UminusCT = Uminus';
        Factor1 = inv(UminusCT * Uminus);
        Factor2 = UminusCT * Uplus;
        A = Factor1 * Factor2;
        e = eig(A);
        for j = 1:Rank
            Zeta(j) = (log(e(j)) * (N));
        end
    end