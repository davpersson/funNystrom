function plotter_low_rank(filename)

%Function to produce Figure 2

addpath('Results')
load(filename)

figure(1)
semilogy(rank_list,trace_optimal/trfA,'k-o','LineWidth',4,'MarkerSize',10)
hold on
semilogy(rank_list,max(mean(trace_error(:,:,1),2)/trfA,1e-16),'b-*','LineWidth',3,'MarkerSize',10)
semilogy(rank_list,max(mean(trace_error(:,:,2),2)/trfA,1e-16),'b-^','LineWidth',3,'MarkerSize',10)
semilogy(rank_list,max(mean(trace_error_exact(:,:,1),2)/trfA,1e-16),'r-x','LineWidth',2,'MarkerSize',10)
semilogy(rank_list,max(mean(trace_error_exact(:,:,2),2)/trfA,1e-16),'r-s','LineWidth',2,'MarkerSize',10)

xlabel('Rank','Interpreter','latex');
ylabel('Relative nuclear norm error','Interpreter','latex');
legend({'Optimal','funNystrom ($q=1$)','funNystrom ($q=2$)',...
    'Nystrom ($q=1$)','Nystrom ($q=2$)'},'Interpreter','latex')
set(gca,'FontSize',18)
hold off

print(append(filename,'_trace'),'-depsc')

figure(2)
semilogy(rank_list,frobenius_optimal/frobfA,'k-o','LineWidth',4,'MarkerSize',10)
hold on
semilogy(rank_list,max(mean(frobenius_error(:,:,1),2)/frobfA,1e-16),'b-*','LineWidth',3,'MarkerSize',10)
semilogy(rank_list,max(mean(frobenius_error(:,:,2),2)/frobfA,1e-16),'b-^','LineWidth',3,'MarkerSize',10)
semilogy(rank_list,max(mean(frobenius_error_exact(:,:,1),2)/frobfA,1e-16),'r-x','LineWidth',2,'MarkerSize',10)
semilogy(rank_list,max(mean(frobenius_error_exact(:,:,2),2)/frobfA,1e-16),'r-s','LineWidth',2,'MarkerSize',10)

xlabel('Rank','Interpreter','latex');
ylabel('Relative Frobenius norm error','Interpreter','latex');
legend({'Optimal','funNystrom ($q=1$)','funNystrom ($q=2$)',...
    'Nystrom ($q=1$)','Nystrom ($q=2$)'},'Interpreter','latex')
set(gca,'FontSize',18)
hold off

print(append(filename,'_frobenius'),'-depsc')

end