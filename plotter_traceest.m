function plotter_traceest(filename)

%Function to produce Figure 4

addpath('Results')
load(filename)

semilogy(m_list,abs((nystrom_approx-tracefA)/tracefA),'b-*','LineWidth',3,'MarkerSize',10)
hold on
semilogy(m_list,abs((rSVD-tracefA)/tracefA),'r-o','LineWidth',3,'MarkerSize',10)
ylabel('Relative error','Interpreter','latex');
xlabel('Matrix-vector products','Interpreter','latex');
legend({'funNystrom','Subspace iteration'},'Interpreter','latex','Location','best')
set(gca,'FontSize',18)
hold off

print(filename,'-depsc')

end