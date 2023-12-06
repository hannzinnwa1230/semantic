clc;
clear;
close all;

fs = 20e6; % Hz
% pathDelays = [0 1 5 50 100]*1e-9; % sec
% avgPathGains = [0 -5.9 -10 -12.3 -15.6]; % dB
% fD = 2; % Hz
% rayleighchan = comm.RayleighChannel('SampleRate',fs, ...
%     'PathDelays',pathDelays, ...
%     'AveragePathGains',avgPathGains, ...
%     'MaximumDopplerShift',fD,...
%     'Visualization','Impulse and frequency responses');
x =[zeros(1000,1);1;zeros(1000,1)];
pathDelays = [0 1 5 50 100]*1e-9; % sec
avgPathGains = [0 -5.9 -10 -12.3 -15.6]; % dB
fD = 2; % Hz
rayleighchan = comm.RayleighChannel('SampleRate',fs, ...
    'PathDelays',pathDelays, ...
    'AveragePathGains',avgPathGains, ...
    'MaximumDopplerShift',fD);
% tgnChan = wlanTGnChannel('SampleRate',fs, ...
%     'LargeScaleFadingEffect','shadowing');
% y = tgnChan(x);
y = rayleighchan(x);
Y =fftshift(fft(y));
plot(abs(Y));
