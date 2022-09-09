%This script reproduces the figures in Figure 2

clc
clear

rng(0)
addpath('Other')
addpath('results')

%Number of repetitions of experiment
repeats = 1;

%Run the experiment for all specified matrices
tic
for matrix = 1:4
        
    if matrix == 1
        
        %Figure 2(a)
        
        %List of ranks of low-rank approximations
        rank_list = 10:10:200;
        
        %Determine the matrix function
        fscalar = @(x) sqrt(x);
        
        %Specify A and f(A)
        matrix_size = 5000;
        parameter = 3;
        U = gallery('orthog',matrix_size);
        Lambda = diag((1:matrix_size).^(-parameter));
        fLambda = diag((1:matrix_size).^(-parameter/2));
        
        %Specify a filename
        filename = 'results/algebraic3_sqrtm_lowrank';
        
    elseif matrix == 2
        
        %Figure 2(b)
        
        %List of ranks of low-rank approximations
        rank_list = 10:10:200;
        
        %Determine the matrix function
        fscalar = @(x) x./(x+0.01);
        
        %Specify A and f(A)
        matrix_size = 5000;
        parameter = 10;
        U = gallery('orthog',matrix_size);
        Lambda = diag(exp(-(1:matrix_size)/parameter));
        fLambda = diag(fscalar(exp(-(1:matrix_size)/(parameter))));
        
        %Specify a filename
        filename = 'results/exponential10_sd001_lowrank';
        
    elseif matrix == 3
        
        %Figure 2(c)
       
        %List of ranks of low-rank approximations
        rank_list = 10:10:200;
        
        %Determine the matrix function
        fscalar = @(x) sqrt(x);
       
        %Specify A and f(A)
        matrix_size = 5000;
        parameter_cell = {3/2,1};
        A = kernel(randn(1,matrix_size),'matern',parameter_cell);
        [U,Lambda,~] = svd(A,'econ');
        fLambda = diag(fscalar(diag(Lambda)));
        
        %Specify a filename
        filename = 'results/matern32_sqrtm_lowrank';

    elseif matrix == 4
        
        %Figure 2(d)
        
        %List of ranks of low-rank approximations
        rank_list = 10:10:200;
        
        %Determine the matrix function
        fscalar = @(x) x./(x+0.01);
       
        %Specify A and f(A)
        matrix_size = 5000;
        parameter_cell = {5/2,1};
        A = kernel(randn(1,matrix_size),'matern',parameter_cell);
        [U,Lambda,~] = svd(A,'econ');
        fLambda = diag(fscalar(diag(Lambda)));
        
        %Specify a filename
        filename = 'results/matern52_sd001_lowrank';
        
    end
    toc
    
    matrix
    
    %Run test
    low_rank_test(U,Lambda,fscalar,repeats,rank_list,filename)
    
    %Plot the results
    plotter_low_rank(filename);
    
end