%This script reproduces Figure 4

clc
clear

rng(0)
addpath('Other')
addpath('results')

for matrix = 1:2

    if matrix == 1
        
        %Figure 4(a)
        
        %Set up matrix
        matrix_size = 5000;
        A = 100*sparse(diag(0.9.^(1:matrix_size)));
        f = @(X) logm(eye(size(X)) + X);
        fscalar = @(x) log(1+x);
        tracefA = sum(fscalar(diag(A)));
        filename = 'results/exponential_traceest';
        
    elseif matrix == 2
        
        %Figure 4(b)
        
        %Set up matrix
        matrix_size = 5000;
        A = 100*sparse(diag((1:matrix_size).^(-3)));
        f = @(X) logm(eye(size(X)) + X);
        fscalar = @(x) log(1+x);
        tracefA = sum(fscalar(diag(A)));
        filename = 'results/algebraic_traceest';
        
    end

%Set up experiment
m_list = 10:10:200;

%Allocate space for results
rSVD = zeros(1,length(m_list));
nystrom_approx = zeros(1,length(m_list));

%Run experiments
iteration = 0;
for m = m_list
    
    iteration = iteration + 1;
    
    %Nystrom
    [~,S] = nystrom(matrix_size,@(X) A*X, m, 1);
    nystrom_approx(iteration) = sum(fscalar(diag(S)));
    
    %rSVD
    [Q,~] = qr(A*randn(matrix_size,m/2),0);
    rSVD(iteration) = trace(f(Q'*A*Q));
    
end

save(filename,'rSVD','nystrom_approx','m_list','tracefA')

%Plot results
plotter_traceest(filename)
end