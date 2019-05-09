[corrupted_signal, K, N, total_iterations] = config();
%%QE methods
rho = 2.25;
QE_Gamma =0.98;
tensor_iterations = total_iterations + 1200;
[X_nothing, y_nothing, Norm_nothing, iterations_nothing, converged_nothing] = QEnothingADMM(corrupted_signal, K, rho, total_iterations, QE_Gamma);
[X_tensor, y_tensor, Norm_tensor, iterations_tensor, converged_tensor] = QEtensorADMM(corrupted_signal, K, rho, tensor_iterations, QE_Gamma);
[X_all, y_all, Norm_all, iterations_all, converged_all] = QEallADMM(corrupted_signal, K, rho, total_iterations, QE_Gamma);

iter = total_iterations/5;
numLambda=500;
Hf = Hank(corrupted_signal,N);
right = norm(Hf);

lambdaNothing = getLambdaNothing(right,corrupted_signal,K,rho,iter,numLambda);
lambdaTensor = getLambdaTensor(right,corrupted_signal,K,rho,iter,numLambda);
lambdaAll = getLambdaAll(right,corrupted_signal,K,rho,iter,numLambda);

[X_nothingNN, y_nothingNN, Norm_nothingNN, iterations_nothingNN, converged_nothingNN ] = NNnothingADMM(corrupted_signal, rho, total_iterations, lambdaNothing);
[X_tensorNN, y_tensorNN, Norm_tensorNN, iterations_tensorNN, converged_tensorNN ] = NNtensorADMM(corrupted_signal, rho, tensor_iterations, lambdaTensor);
[X_allNN, y_allNN, Norm_allNN, iterations_allNN, converged_allNN ] = NNallADMM(corrupted_signal, rho, total_iterations, lambdaAll);

[X_nothingHC, y_nothingHC, Norm_nothingHC, iterations_nothingHC, converged_nothingHC ] = HCnothingADMM(corrupted_signal, K, rho, total_iterations);
[X_tensorHC, y_tensorHC, Norm_tensorHC, iterations_tensorHC, converged_tensorHC ] = HCtensorADMM(corrupted_signal, K, rho, tensor_iterations);
[X_allHC, y_allHC, Norm_allHC, iterations_allHC, converged_allHC ] = HCallADMM(corrupted_signal, K, rho, total_iterations);
