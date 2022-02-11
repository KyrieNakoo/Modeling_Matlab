% Load data 
load('cortisol.mat'); 

T = 0.5;
fs = 1/T;
% time = linspace(0, length(cortisol)/fs, length(cortisol));
time = (0:length(cortisol)-1)'/fs;
plot(time,cortisol);
xlabel('time (hr)'); 
ylabel('Plasma cortisol');

%% Correlation
[c,lag] = xcorr(cortisol - mean(cortisol), 'coeff');
figure; plot(lag/fs, c);
xlabel('Lag (hr)'); ylabel('Auto-correlation');

%% DFT
T = 0.5;
fs = 1/T;
% time = linspace(0, length(cortisol)/fs, length(cortisol));
time = (0:length(cortisol)-1)'/fs;
nfft = length(cortisol);

y = cortisol-mean(cortisol);
Y = fft(y, nfft);
Y = Y(1:nfft/2+1);
Y = (2/length(y)) * Y;
Ymag = abs(Y);
Ymag(1) = Ymag(1)/2;

f = linspace(0,fs/2,nfft/2+1)';
plot(f, Ymag);
xlabel('Frequency (Hz)'); ylabel('|Cortisol (f)|');



