function [error_lanczos,matvecs] = accuracy_test_lanczos(U,Lambda,fLambda,fscalar,r,q,error_measure)

%SECTION 1: Set parameters
fprintf('Starting Nystrom with exact matvecs \n')
tolerance = 1.1;
n_it_max = 50; %Maximum number of lanczos iterations
step = 5;

%SECTION 2: Compute the low rank approximation using exact matvecs
matrix_size = size(U,1);
Omega = randn(matrix_size,r);
fAfun = @(X) U*(fLambda*(U'*X));
Afun = @(X) U*(Lambda*(U'*X));

[Uhat,Lambdahat] = nystrom(matrix_size,fAfun,r,q,Omega);

%Compute error
error_exact = error_measure_fun(U*fLambda*U',Uhat*Lambdahat*Uhat',error_measure)

%SECTION 3: Compute low rank approximation with approximate matrix-vector
%products
error_lanczos = inf;
n_it = 0;

fprintf('Starting testing how many Lanczos iterations we need \n')

flag = 1;

while (n_it < n_it_max)&&(error_lanczos > tolerance*error_exact)&&(flag)
    
    error_lanczos_old = error_lanczos;
    
    n_it = n_it + step;
    fprintf('Lanczos iterations %i \n',n_it)
    matvecs = n_it*q*r;
    
    fAfun_block_lanczos = @(X) block_lanczos(Afun,X,fscalar,n_it);
    [Uhat,Lambdahat] = nystrom(matrix_size,fAfun_block_lanczos,r,q,Omega);
    
    %Compute error
    error_lanczos = error_measure_fun(U*fLambda*U',Uhat*Lambdahat*Uhat',error_measure)
    
    if (n_it > 1) && ((abs(error_lanczos-error_lanczos_old)/error_lanczos_old) < 0.01)
        
        flag = 0;
        
    end
    
end
    
end