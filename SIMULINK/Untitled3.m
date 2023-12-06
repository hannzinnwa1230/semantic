clc;
clear;
close all;

load lena44.mat
load lena46.mat
load lena48.mat
load peppers44.mat
load peppers46.mat
load peppers48.mat
load baboon44.mat
load baboon46.mat
load baboon48.mat
load airplane44.mat
load airplane46.mat
load airplane48.mat

figure(1)
plot(1:50,lena44(end,:),"--r");
hold on;
plot(1:50,lena46(end,:),"-r");
hold on; 
plot(1:50,lena48(end,:),"-.r");
legend('γ=44','γ=46','γ=48');
grid on
xlabel('遗传代数')
ylabel('最优适应度')


figure(2)
plot(1:50,peppers44(end,:),"--r");
hold on;
plot(1:50,peppers46(end,:),"-r");
hold on; 
plot(1:50,peppers48(end,:),"-.r");
legend('γ=44','γ=46','γ=48');
grid on
xlabel('遗传代数')
ylabel('最优适应度')


figure(3)
plot(1:50,baboon44(end,:),"--r");
hold on;
plot(1:50,baboon46(end,:),"-r");
hold on; 
plot(1:50,baboon48(end,:),"-.r");
legend('γ=44','γ=46','γ=48');
grid on
xlabel('遗传代数')
ylabel('最优适应度')


figure(4)
plot(1:50,airplane44(end,:),"--r");
hold on;
plot(1:50,airplane46(end,:),"-r");
hold on; 
plot(1:50,airplane48(end,:),"-.r");
legend('γ=44','γ=46','γ=48');
grid on
xlabel('遗传代数')
ylabel('最优适应度')




