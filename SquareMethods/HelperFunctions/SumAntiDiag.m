function out=SumAntiDiag(A,N)
  out=zeros(1,2*N+1);for j=1:2*N+1, out(j)=sum(diag(flipud(A),-N-1+j));
end
    