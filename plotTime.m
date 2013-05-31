clc, clear all

% plot timing.txt data
filename = 'stout/timing_830268.txt'
fid = fopen(filename);
data = textscan(fid, '%d %f %f %f %f %f %f');
fclose(fid);

% T-ref
T1=data{6}(1)*8; % T1 is the execution time of the sequential algorithm
% Speedup
speedup = T1./data{6}; % T1 divided by the execution time of the parallel algorithm with p processors

%%
figure(1), clf
% set(gca,'XScale','log','YScale','lin')

hold on
plot(data{1},speedup,'x','MarkerSize',12) % obeserved
plot(data{1},data{1},'-k','LineWidth',2) % ideal 
% plot(data{1},TD,'-xk')
hold off
grid on
legend(sprintf('nelX x nelY x nelZ, steps: %d x %d x %d, %d',data{2}(1),data{3}(1),data{4}(1),data{5}(1)),'ideal','location','nw')
title('Scalability of 2D implementation')
xlabel('number of processors')
ylabel('speedup, T_1/T')
ylim([0 max(speedup)])

print(gcf,'-dpng','timePlot.png')