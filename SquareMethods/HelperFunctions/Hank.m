function out=Hank(a,N)
    out=hankel(a(1:N+1),a(N+1:2*N+1));
end