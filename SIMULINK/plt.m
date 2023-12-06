clc;
clear;
close all;

load test44.mat

plot(1:MAXGEN,trace(end,:),"--r");
hold on;

load test46.mat

plot(1:MAXGEN,trace(end,:),"-r");
hold on;

load test48.mat

plot(1:MAXGEN,trace(end,:),"-.r");
hold on;

legend('γ=44','γ=46','γ=48');
grid on
xlabel('遗传代数')
ylabel('最优适应度')

% xlabel('遗传代数')
% ylabel('性能变化')
% title('进化过程')
% bestY = trace(end,end);
% bestX=trace(1:end-1,end);
% fprintf(['最优嵌入因子:\nX=',num2str(bestX'),'\n最优性能err=',num2str(bestY),'\n'])

% figure(1)
% plot(1:50,lena44(end,:),"--r");
% hold on;
% plot(1:50,lena46(end,:),"-r");
% hold on; 
% plot(1:50,lena48(end,:),"-.r");
% legend('γ=44','γ=46','γ=48');
% grid on
% xlabel('遗传代数')
% ylabel('最优适应度')