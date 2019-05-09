%%%
% This file runs the ADMM algorithm using the square case of the Nuclear Norm without any tensors or misfit weights as discussed in paper
%
% {params} 
% f - signal vector input of length 2 * N - 1
% K - number of complex exponentials used
% rho - value of rho, can affect rate of convergence in ADMM
% iterations- number of iterations to run before breaking ADMM
%
% {returns}
% matrix X of shape (2*N+1) * (2*N+1),
% vector y of length 2 * N - 1,
% differenceNorm - norm of (X - H(y)), where H(y) is the Hankel matrix created from the vector y,
% iterationsTaken - number of iterations taken to converge,
% converged - boolean whether method converged
%
%%% 
function [X,y,differenceNorm,iterationsTaken,converged] = HCnothingADMM(f,K,rho,iterations)
    N=(length(f)-1)/2; %f has length 2N+1
    y=zeros(1, 2*N+1);
    X=ones(N+1, N+1);
    L=zeros(N+1, N+1);
    iterationsTaken = 0;
    
    u = ones(N+1,1);
    v = ones(1,N+1);
    W=u*v;  
    t=SumAntiDiag(W,N);
    
    while norm(X - Hank(y, N)) > 10^-12 && iterationsTaken<iterations
        cutvector=zeros(N + 1, 1);
        cutvector(1:K)=1;
        [U, S, V] = svd((Hank(y, N) - L));
        s = diag(S);
        s = s.* cutvector;
        X = U * diag(s) * V';
        % y update step
        y = (t.* f + SumAntiDiag(rho * (X + L), N))./(t + rho*t);
        % Lambda update Step
        L = L + (X - Hank(y, N));
        iterationsTaken = iterationsTaken + 1;
    end

    differenceNorm = norm(X - Hank(y, N));
    converged = (differenceNorm > 10^-12);
end