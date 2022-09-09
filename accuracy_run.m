%This script reproduces the figures in Figure 1

clc
clear

rng(0)

addpath('Other')
addpath('results')

%Number of subspace iterations
q = 1;

%Run the experiment for all specified matrices
for matrix = 1:4
        
    if matrix == 1
        
        %Figure 1(a)
        
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
        
        %Specify how we measure the error
        error_measure = 'frobenius_norm';
        
        %Specify a filename
        filename = 'results/algebraic3_sqrtm_accuracy';
        
    elseif matrix == 2
        
        %Figure 1(b)
        
        %List of ranks of low-rank approximations
        rank_list = 10:10:200;
        
        %Determine the matrix function
        fscalar = @(x) x./(x+1);
        
        %Specify A and f(A)
        matrix_size = 5000;
        parameter = 10;
        U = gallery('orthog',matrix_size);
        Lambda = 10*diag(exp(-(1:matrix_size)/parameter));
        fLambda = diag(fscalar(diag(Lambda)));
        
        %Specify how we measure the error
        error_measure = 'nuclear_norm';
        
        %Specify a filename
        filename = 'results/exponential10_sd1_accuracy'
        
    elseif matrix == 3
        
        %Figure 1(c)
        
        %List of ranks of low-rank approximations
        rank_list = 10:10:50;
        
        %Determine the matrix function
        fscalar = @(x) log(1+x);
       
        %Specify A and f(A)
        matrix_size = 5000;
        parameter_cell = {0.1};
        A = kernel(randn(1,matrix_size),'square_exponential',parameter_cell);
        [U,Lambda,~] = svd(A,'econ');
        fLambda = diag(fscalar(diag(Lambda)));
        
        ;%Specify how we measure the error
        error_measure = 'frobenius_norm';
        
        %Specify a filename
        filename = 'results/SE01_logm_accuracy';

    elseif matrix == 4
        
        %Figure 1(d)
        
        %List of ranks of low-rank approximations
        rank_list = 10:10:100;
        
        %Determine the matrix function
        fscalar = @(x) log(x+1);
        
        %Specify A
        kappa = 0.01;
        sigma = 1;
        lambda = 1;
        T_list = [1,1.5,2];
        A = uq_laplaceeig(kappa,sigma,lambda,T_list);
        [U,Lambda,~] = svd(A,'econ');
        fLambda = diag(fscalar(diag(Lambda)));
        matrix_size = size(A,1);
        
        %Specify how we measure the error
        error_measure = 'nuclear_norm';
        
        %Specify a filename
        filename = 'results/pde_lambda=1_logm_accuracy';
        
    end
    
    %Allocate space for results
    error_funnystrom_list = zeros(length(rank_list),1);
    error_lanczos_list = zeros(length(rank_list),1);
    matvec_funnystrom_list = zeros(length(rank_list),1);
    matvec_lanczos_list = zeros(length(rank_list),1);

    fprintf('Testing on matrix %i \n',matrix)
    iteration = 0;
    for r = rank_list

        fprintf('Rank %i \n',r)
        iteration = iteration + 1;
        
        %Run test
        [error_funnystrom,error_lanczos,matvec_funnystrom,matvec_lanczos] ...
            = accuracy_test(U,Lambda,fLambda,fscalar,r,q,error_measure);

        error_funnystrom_list(iteration) = error_funnystrom;
        error_lanczos_list(iteration) = error_lanczos;
        matvec_funnystrom_list(iteration) = matvec_funnystrom;
        matvec_lanczos_list(iteration) = matvec_lanczos;
        
    end
    
    %Save results
    save(filename,'error_funnystrom_list','error_lanczos_list',...
        'matvec_funnystrom_list','matvec_lanczos_list',...
        'fLambda','error_measure')
    
    %Plot the results
    plotter_accuracy(filename);
    
end