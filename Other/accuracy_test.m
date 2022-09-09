function [error_funnystrom,error_lanczos,...
    matvec_funnystrom,matvec_lanczos] = accuracy_test(U,Lambda,fLambda,fscalar,r,q,error_measure)

%SECTION 1: Set up experiment
%Compute important quantities
fprintf('Starting funNystrom \n')
matrix_size = size(U,1);
Afun = @(X) U*Lambda*U'*X;

%SECTION 2: Compute errors for funNystrom
matvec_funnystrom = q*r;

[Uhat,Lambdahat] = nystrom(matrix_size,Afun,r,q); 
fLambdahat = diag(fscalar(diag(Lambdahat)));

%Compute the error
error_funnystrom = error_measure_fun(U*fLambda*U',Uhat*fLambdahat*Uhat',error_measure)

fprintf('Starting Nystrom \n')

%SECTION 3: Compute errors for Nystrom with Lanczos
[error_lanczos,matvec_lanczos] = accuracy_test_lanczos(U,Lambda,fLambda,fscalar,r,q,error_measure);

end