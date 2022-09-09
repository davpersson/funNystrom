function t = nystrompp(matrix_size,Afun,m1,m2)

%Implementation of Nystrom++

%Low rank approximation phase
[U,S] = nystrom(matrix_size,Afun,m1,1);
t1 = trace(S);

%Stochastic trace estimation phase
Omega = randn(matrix_size,m2);
UtOmega = U'*Omega;
t2 = (trace(Omega'*Afun(Omega))-trace(UtOmega'*S*UtOmega))/m2;

%Return result
t = t1 + t2;

end