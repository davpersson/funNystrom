function plotter_time(filename)

%Function to produce Figure 3

addpath('Results')
load(filename)

figure('Renderer', 'painters', 'Position', [10 10 900 300])
plot(N_list,time_lanczos./time_lowrank,'k-*','LineWidth',3,'MarkerSize',10)
hold on
ylabel('Speed up factor','Interpreter','latex');
xlabel('$N$','Interpreter','latex');
set(gca,'FontSize',18)
hold off

print(filename,'-depsc')

end