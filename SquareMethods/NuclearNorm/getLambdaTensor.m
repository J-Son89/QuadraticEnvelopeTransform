
function lambda=getLambdaTensor(right,corrupted_signal,K,rho,iter,numLambda)
lamba_iter = 0 ;
left=0;
for i=1:numLambda
     mid = (left+right)/2;
    [XX,y,~,~,~]=NNtensorADMM(corrupted_signal,rho,floor(iter),mid);
     if rank(XX)>K %%rank too high decrease lambda
         left =mid;
     elseif rank(XX)<K %%rank too high decrease lambda
         right =mid;
     else %% comes in here once rank is correct
         lambda = mid;
         right=mid; %% sets right endpoint to current value of lambda where rank is correct
         lamba_iter = i;%% number of iterations to get oldlambda
         break;        %% misfit increases as lambda decreases. want right to be as small as possible
     end 
 end
 current_misfit = norm(corrupted_signal - y); 
 i=0;
 while i<50 
     mid = (left+right)/2; %% uncertain whether using this mid point value for lambda returns rank K
     [XX,y,~,~,~]=NNtensorADMM(corrupted_signal,rho,floor(iter),mid); 
     if rank(XX)>K %%rank from using midpoint too high => move left end point up
        left =mid;
     else
        right = mid; %%otherwise rank using mid point is K so we can move our right end point lower
        lambda =right; %% set lambda to right end point, we can be sure method returns correct rank
     end
     
     if current_misfit -norm(corrupted_signal - y) < 0.1
         break
     else 
          current_misfit = norm(corrupted_signal - y);
     end
     i=i+1;
 end
    lamba_iter = lamba_iter +i; %% number of iterations to get lambda
end