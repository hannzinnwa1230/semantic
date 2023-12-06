clc;
clear;
close all;

% load MCS1.mat
% PER1 = PER;
% SNRrecord1 = SNRrecord;
% 
% load MCS4.mat
% PER2 = PER;
% SNRrecord2 = SNRrecord;
% 
% 
% x = [8:0.2:16];
% y1 = ones(1,41);
% y2 = zeros(1,41);
% y1(end-20:end) = PER1;
% y2(1:23) = PER2;
% plot(x,y1,'--')
% hold on
% plot(x,y2,'-')
% legend('64QAM 1/2码率','16QAM 1/2码率')
% ylabel('误包率')
% xlabel('信噪比/dB')

load MCS1.mat
PER1 = PER;
SNRrecord1 = SNRrecord;

load MCS6.mat
PER2 = PER;
SNRrecord2 = SNRrecord;

x = [12:0.2:19];
y2 = ones(1,36);
y1 = zeros(1,36);
y1(1:21) = PER1;
y2(end-25:end) = PER2;
plot(x,y1,'--')
hold on
plot(x,y2,'-')
legend('64QAM 1/2码率','64QAM 2/3码率')
ylabel('误包率')
xlabel('信噪比/dB')