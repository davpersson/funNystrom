function [U,S] = nystrom(matrix_size,Afun,r,q,varargin)

%Implementation of the Nystrom approximation

%Generate sketch matrix
if nargin < 5
    
    [Omega,~] = qr(randn(matrix_size,r),0);
    
elseif nargin == 5
    
    Omega = varargin{1};
    [Omega,~] = qr(Omega,0);
    
else
    
    error('Too many input arguments')
    
end

iteration = 0;

%Subspace iteration
while iteration < q-1
    
    [Omega,~] = qr(Afun(Omega),0);
    iteration = iteration + 1;
    
end

%Compute matrix products with A
Y = Afun(Omega);

%Regularizaton
[V,D,~] = svd(Omega'*Y,'econ');
D(D < 5e-16*D(1,1)) = 0;
B = Y*(V*pinv(diag(sqrt(diag(D))))*V');
[U,Shat,~] = svd(B,'econ');
S = Shat^2;

end