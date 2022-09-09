function plotter_accuracy(filename)

%Figure to produce Figure 1

addpath('Results')
load(filename)

if strcmp(error_measure,'nuclear_norm')
    
    normfA = trace(fLambda);
    x_label = 'Relative nuclear norm error';
    
elseif strcmp(error_measure,'frobenius_norm')
    
    normfA = norm(diag(fLambda));
    x_label = 'Relative Frobenius norm error';
    
else
    
    error('error_measure must be nuclear_norm or frobenius_norm')
    
end

loglog(error_funnystrom_list/normfA,matvec_funnystrom_list,'b-*','LineWidth',3,'MarkerSize',10)
hold on
loglog(error_lanczos_list/normfA,matvec_lanczos_list,'r-x','LineWidth',3,'MarkerSize',10)

ylabel('Matrix-vector products','Interpreter','latex');
xlabel(x_label,'Interpreter','latex');
legend({'funNystrom','Nystrom approximation'},'Interpreter','latex','Location','best')
set(gca,'FontSize',18)
hold off

print(filename,'-depsc')

end