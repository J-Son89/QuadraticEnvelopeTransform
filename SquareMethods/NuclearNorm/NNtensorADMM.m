%%%
% This file runs the ADMM algorithm using the square case of the Nuclear Norm without any tensors or misfit weights as discussed in paper
%
% {params} 
% f - signal vector input of length 2 * N - 1
% rho - value of rho, can affect rate of convergence in ADMM
% iterations- number of iterations to run before breaking ADMM
% lambda - parameter for nuclear norm method
%
% {returns}
% matrix X of shape (2*N+1) * (2*N+1),
% vector y of length 2 * N - 1,
% differenceNorm - norm of (X - H(y)), where H(y) is the Hankel matrix created from the vector y,
% iterationsTaken - number of iterations taken to converge,
% converged - boolean whether method converged
%
%%% 
function [X,y,differenceNorm,iterationsTaken,converged] = NNtensorADMM(f,rho,iterations, lambda)
    N=(length(f)-1)/2; %f has length 2N+1
    y=zeros(1, 2*N+1);
    X=ones(N+1, N+1);
    L=zeros(N+1, N+1);
    u = zeros(N+1,1);
    v = zeros(1,N+1);
    for j = 1:N+1
        u(j) = 1/ sqrt( (N/2) +1 - abs(j-1-N/2));
        v(j) = 1/ sqrt( (N/2)+1 - abs(j-1-N/2));
    end
    W=u*v;  
    t=SumAntiDiag(W,N);
    iterationsTaken = 0;

    while norm(X - Hank(y, N)) > 10^-12 && iterationsTaken<iterations
        [U,S,V] = svd((W.^(1/2)).*(Hank(y, N) - L));
        s = diag(S);
        s = max(s - lambda, 0);
        %Xupdate step
        X = U * diag(s) * V';
        X = X./(W.^(1/2));
        % y update step
        y = (t.* f + SumAntiDiag(rho * W.*(X + L), N))./(t + rho * t);
        % Lambda update Step
        L = L + (X - Hank(y, N));
        iterationsTaken = iterationsTaken + 1;
    end

    differenceNorm = norm(X - Hank(y, N));
    converged = (differenceNorm > 10^-12);
end