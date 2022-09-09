function t = funnystrompp(matrix_size,Afun,fAfun,r,m,q,fscalar)

%Implementation of Algorithm 2

%Low rank approximation phas
[U,S] = nystrom(matrix_size,Afun,r,q);
t1 = sum(fscalar(diag(S)));

%Stochastic trace estimation phase
Omega = randn(matrix_size,m); Y = U'*Omega;
t2 = (trace(Omega'*fAfun(Omega)) - trace(Y'*diag(fscalar(diag(S)))*Y))/m;

%Return result
t = t1 + t2;

end