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
        
        d = diag(Lambda); fd = diag(fLambda);
        F = @(l,k,q) sqrt((1 + ((d(k+1)/d(k))^(2*(q-1)))*(k/(l-k-1))))*norm(fd((k+1):end));
        normoption = 'fro';
        
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
        
        d = diag(Lambda); fd = diag(fLambda);
        F = @(l,k,q) (1 + ((d(k+1)/d(k))^(2*(q-1)))*(k/(l-k-1)))*sum(fd((k+1):end));
        normoption = 'trace';
        
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
        
        d = diag(Lambda); fd = diag(fLambda);
        F = @(l,k,q) sqrt((1 + ((d(k+1)/d(k))^(2*(q-1)))*(k/(l-k-1))))*norm(fd((k+1):end));
        normoption = 'fro';

        
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
        
        d = diag(Lambda); fd = diag(fLambda);
        F = @(l,k,q) (1 + ((d(k+1)/d(k))^(2*(q-1)))*(k/(l-k-1)))*sum(fd((k+1):end));
        normoption = 'trace';
        
    end
    
    bound = zeros(length(rank_list),2);
    iteration = 0;
    
    for r = rank_list
        
        r
        
        iteration = iteration + 1;
        
        FF = @(k) F(r,k,1);
        bound(iteration,1) = minfinder(FF,r);
        
        FF = @(k) F(r,k,2);
        bound(iteration,2) = minfinder(FF,r);
        
    end
    
    plotter(filename,bound,normoption)
    
    
end


function plotter(filename,bound,normoption)

addpath('Results')
load(filename)

if strcmp(normoption,'trace')
    semilogy(rank_list,trace_optimal/trfA,'k-o','LineWidth',4,'MarkerSize',10)
    hold on
    semilogy(rank_list,max(mean(trace_error(:,:,1),2)/trfA,1e-16),'b-*','LineWidth',3,'MarkerSize',10)
    semilogy(rank_list,max(mean(trace_error(:,:,2),2)/trfA,1e-16),'b-^','LineWidth',3,'MarkerSize',10)
    semilogy(rank_list,max(mean(trace_error_exact(:,:,1),2)/trfA,1e-16),'r-x','LineWidth',2,'MarkerSize',10)
    semilogy(rank_list,max(mean(trace_error_exact(:,:,2),2)/trfA,1e-16),'r-s','LineWidth',2,'MarkerSize',10)
    
    semilogy(rank_list,bound(:,1)/trfA,'b--*','LineWidth',2,'MarkerSize',10)
    semilogy(rank_list,bound(:,2)/trfA,'b--^','LineWidth',2,'MarkerSize',10)

    xlabel('Rank','Interpreter','latex');
    ylabel('Relative nuclear norm error','Interpreter','latex');
    legend({'Optimal','funNystrom ($q=1$)','funNystrom ($q=2$)',...
        'Nystrom ($q=1$)','Nystrom ($q=2$)',...
        'Theoretical bound ($q=1$)','Theoretical bound ($q=2$)'},'Interpreter','latex','Location','best')
    set(gca,'FontSize',18)
    hold off

    print(append(filename,'_trace_bound'),'-depsc')

else
    semilogy(rank_list,frobenius_optimal/frobfA,'k-o','LineWidth',4,'MarkerSize',10)
    hold on
    semilogy(rank_list,max(mean(frobenius_error(:,:,1),2)/frobfA,1e-16),'b-*','LineWidth',3,'MarkerSize',10)
    semilogy(rank_list,max(mean(frobenius_error(:,:,2),2)/frobfA,1e-16),'b-^','LineWidth',3,'MarkerSize',10)
    semilogy(rank_list,max(mean(frobenius_error_exact(:,:,1),2)/frobfA,1e-16),'r-x','LineWidth',2,'MarkerSize',10)
    semilogy(rank_list,max(mean(frobenius_error_exact(:,:,2),2)/frobfA,1e-16),'r-s','LineWidth',2,'MarkerSize',10)

    semilogy(rank_list,bound(:,1)/frobfA,'b--*','LineWidth',2,'MarkerSize',10)
    semilogy(rank_list,bound(:,2)/frobfA,'b--^','LineWidth',2,'MarkerSize',10)
    
    xlabel('Rank','Interpreter','latex');
    ylabel('Relative Frobenius norm error','Interpreter','latex');
    legend({'Optimal','funNystrom ($q=1$)','funNystrom ($q=2$)',...
        'Nystrom ($q=1$)','Nystrom ($q=2$)',...
        'Theoretical bound ($q=1$)','Theoretical bound ($q=2$)'},'Interpreter','latex','Location','best')
    set(gca,'FontSize',18)
    hold off
    
    print(append(filename,'_frobenius_bound'),'-depsc')

end

end

function minimum = minfinder(FF,l)

minimum = Inf;

for i = 1:(l-2)
    
    f = FF(i);
    
    if minimum > f
        
        minimum = f;
        
    end
    
end

end