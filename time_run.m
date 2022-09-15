%This script reproduces the Figure 3

clc
clear

rng(0)
addpath('Other')
addpath('results')

%Create matrix
disp('Create matrix')
matrix_size = 10000;
Q = gallery('orthog',matrix_size,1);
D = diag(exp(-(1:matrix_size)));
Afun = @(X) Q*(diag(D).*(Q'*X));

%Determine the matrix function
f = @(X) sqrtm(X);
fscalar = @(x) sqrt(x);

%Specify a filename
filename = 'results/comptimes';

%Specify parameters
phat = 21;
rhat = 14;
    
%Allocate space to save results
N_list = 100:100:1000;
N_list = 10:10:100;
time_lanczos = zeros(1,length(N_list));
time_lowrank = zeros(1,length(N_list));

%Generate matrix to compute f(A) *  Omega
Omega = eye(matrix_size);

iteration = 0;
for N = N_list
    
    fprintf("%i / %i\n",N,N_list(end))
    
    iteration = iteration + 1;
    func = @() lowrank_matvec(matrix_size, Afun, fscalar, rhat, Omega(:,1:N));
    time_lowrank(iteration) = timeit(func);

    func = @() lanczos_matvec(matrix_size, Afun, fscalar, phat,Omega(:,1:N));
    time_lanczos(iteration) = timeit(func);

end


save(filename,'N_list','time_lanczos','time_lowrank');
%Plot the results
plotter_time(filename);

function lanczos_matvec(matrix_size, Afun, fscalar, p, Omega)
Y = block_lanczos(Afun,Omega,fscalar,p);
end

function lowrank_matvec(matrix_size, Afun, fscalar, rank, Omega)
[U,S] = nystrom(matrix_size, Afun, rank, 1);
Y = U*(fscalar(diag(S)).*(U'*Omega));
end