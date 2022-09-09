function plotter_nystrompp(filename)

%Function to produce Figure 5

addpath('Other')
addpath('results')

load(filename);
matvecs = matvecsA;

relative_errors_fhpp = sort(abs((funNystrompp(:,:,1)-tracefA)/tracefA));
relative_errors_fhpp1 = sort(abs((funNystrompp(:,:,2)-tracefA)/tracefA));
relative_errors_npp = sort(abs((Nystrompp-tracefA)/tracefA));

mean_error_fhpp = mean(abs((funNystrompp(:,:,1)-tracefA)/tracefA));
mean_error_fhpp1 = mean(abs((funNystrompp(:,:,2)-tracefA)/tracefA));
mean_error_npp = mean(abs((Nystrompp-tracefA)/tracefA));

prctile_10_fhpp = relative_errors_fhpp(10,:);
prctile_90_fhpp = relative_errors_fhpp(90,:);
prctile_10_npp = relative_errors_npp(10,:);
prctile_90_npp = relative_errors_npp(90,:);


semilogy(matvecs,prctile_10_npp,'r--')
hold on
semilogy(matvecs,prctile_90_npp,'r--')
patch([matvecs fliplr(matvecs)],[prctile_10_npp...
    fliplr(prctile_90_npp)],'r')

semilogy(matvecs,prctile_10_fhpp,'b--')
semilogy(matvecs,prctile_90_fhpp,'b--')
patch([matvecs fliplr(matvecs)],[prctile_10_fhpp...
    fliplr(prctile_90_fhpp)],'b')


alpha(.1)

h(1)=semilogy(matvecs,mean_error_fhpp1,'k--x','LineWidth',3,'MarkerSize',10);
h(2)=semilogy(matvecs,mean_error_fhpp,'b-*','LineWidth',3,'MarkerSize',10);
h(3)=semilogy(matvecs,mean_error_npp,'r-o','LineWidth',3,'MarkerSize',10);

xlabel('Matrix-vector products with $\mathbf{A}$','interpreter','latex')
ylabel('Relative error','interpreter','latex')
legend(h,{'funNystrom++ (q = 1)','funNystrom++ (q = 2)','Nystrom++'},'interpreter','latex')

set(gca,'FontSize',20)
hold off

print(filename,'-depsc')

end