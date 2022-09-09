%This script reproduces Figure 5

clc
clear

rng(0)
addpath('Other')
addpath('results')

for matrix = 1:2

    if matrix == 1
        
        %Set up matrix
        matrix_size = 5000;
        D = 100*sparse(diag((1:matrix_size).^(-2)));
        Q = gallery('orthog',matrix_size,1);
        f = @(X) logm(eye(size(X)) + X);
        fscalar = @(x) log(1+x);
        tracefA = sum(fscalar(diag(D)));
        
        %Set up experiment
        %Parameters that can be chosen
        lanczos_iterations = 10;
        m_list = 12:12:120; %Matvecs with f(A)
        repeats = 100;
        q = 2;
        
        %Parameters that are determined by the above parameters or should
        %not be changed
        matvecsA = lanczos_iterations*m_list;
        r_list = lanczos_iterations*m_list/(2*q);
        Afun = @(X) Q*(D*(Q'*X));
        fAfun = @(X) block_lanczos(Afun, X, fscalar, lanczos_iterations); 
        
        filename = 'results/algebraic_c=2_logm_npp';
        
    elseif matrix == 2
        
        %Set up matrix
        matrix_size = 5000;
        D = sparse(diag(exp(-(1:matrix_size)/100)));
        %D = sparse(diag((1:matrix_size)).^(-3));
        Q = gallery('orthog',matrix_size,1);
        f = @(X) logm(eye(size(X)) + X);
        fscalar = @(x) x./(x+0.1);
        tracefA = sum(fscalar(diag(D)));
        
        %Set up experiment
        %Parameters that can be chosen
        lanczos_iterations = 10;
        m_list = 12:12:120; %Matvecs with f(A)
        repeats = 100;
        q = 2;
        
        %Parameters that are determined by the above parameters or should
        %not be changed
        matvecsA = lanczos_iterations*m_list;
        r_list = lanczos_iterations*m_list/(2*q);
        Afun = @(X) Q*(D*(Q'*X));
        fAfun = @(X) block_lanczos(Afun, X, fscalar, lanczos_iterations); 
        
        filename = 'results/exponential_s=100_sd01_npp';
        
    end

    %Allocate space for results
    funNystrompp = zeros(repeats,length(m_list),2);
    Nystrompp = zeros(repeats,length(m_list));

    %Run experiments
    for iteration = 1:length(m_list)

        m = m_list(iteration);
        r = r_list(iteration);
        iteration/length(m_list)

        for repetition = 1:repeats

            Nystrompp(repetition,iteration) = nystrompp(matrix_size,fAfun,m/2,m/2);
            funNystrompp(repetition,iteration,1) = funnystrompp(matrix_size,Afun,fAfun,r,m/2,q,fscalar);
            funNystrompp(repetition,iteration,2) = funnystrompp(matrix_size,Afun,fAfun,q*r,m/2,1,fscalar);

        end

    end
    
    save(filename,'funNystrompp','Nystrompp','tracefA','matvecsA')
    
    plotter_nystrompp(filename)
    
end