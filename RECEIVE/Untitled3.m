clear all
load txWaveform
rxWaveform = resample1(txWaveform,fs,fs*osf); % 降采样
rxWaveform = resample(txWaveform,fs,fs*osf); % 降采样