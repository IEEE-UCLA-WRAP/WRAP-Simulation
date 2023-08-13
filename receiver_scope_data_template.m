%% Oscilloscope Data Template
%% Constants
% to do on MC: calibrate to find the sampling frequency
%Sampling frequency
Fs = 5000000;

%Number of samples
N=50000;

%Time for samples
Tmax = N/Fs;

% from other sim, actual frequency ~9.807e5
%% Reading Scope Data
[x,y] = importAgilentBin("scope_4_1.bin");
x = transpose(x);
y = transpose(y);

% Make time start at 0
x = x+Tmax/2;
figure;
plot (x,y);
legend("Scope Signal");

%% Implement your receiver code below





