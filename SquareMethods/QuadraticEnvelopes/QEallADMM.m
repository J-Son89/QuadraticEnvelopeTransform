%%%
% This file runs the ADMM algorithm using the square case of the Quadratic Envelope with tensors and misfit weights as discussed in paper
%
% {params} 
% f - signal vector input of length 2 * N - 1
% K - number of complex exponentials used
% rho - value of rho, can affect rate of convergence in ADMM
% iterations- number of iterations to run before breaking ADMM
% gama (gamma) - alters the convexity - < 1 is convex
%
% {returns}
% matrix X of shape (2*N+1) * (2*N+1),
% vector y of length 2 * N - 1,
% differenceNorm - norm of (X - H(y)), where H(y) is the Hankel matrix created from the vector y,
% iterationsTaken - number of iterations taken to converge,
% converged - boolean whether method converged,
%
%%% 
function [X,y,differenceNorm,iterationsTaken,converged] = QEallADMM(f,K,rho,iterations,gama)
    N=(length(f)-1)/2; %f has length 2N+1
    
    u = zeros(N+1,1);
    v = zeros(1,N+1);
    for j = 1:N+1
        u(j) = 1/ sqrt( (N/2) +1 - abs(j-1-N/2));
        v(j) = 1/ sqrt( (N/2)+1 - abs(j-1-N/2));
    end
    W=u*v;
    t=SumAntiDiag(W,N);
    tUseMisfitW=t;
    
    y=zeros(1, 2*N+1);
    X=ones(N+1,N+1);
    L=zeros(N+1,N+1);
    
    iterationsTaken=0;
    
    while norm(X-Hank(y, N)) > 10^-12 && iterationsTaken < iterations
        [U, S, V] = svd(( W.^(1/2)).*( Hank(y, N) - L ));
        s = diag(S);
        s = Prox2(K, gama,rho,s');
        X = U * diag(s) * V';
        X = X./(W.^( 1/2 ));
        % y update step
        y = (t.* f./tUseMisfitW + SumAntiDiag(rho * W.*(X + L), N))./(t./tUseMisfitW + rho * t);
        %Lambda update step
        L = L + (X - Hank(y, N));
        iterationsTaken = iterationsTaken + 1;
        
    end
    differenceNorm = norm(X - Hank(y, N));
    converged = (differenceNorm > 10^-12);
end