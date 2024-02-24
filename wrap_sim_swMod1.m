close all;
clear;
clc;

%% Setup

% Values that we're giving you
Fs = 4e6; % Sampling rate (samples/sec)
Rs = 50e3; % Symbol rate in symbols/sec (baud)
Tmax = 0.1; % Max time for the simulation (sec)
fc = 1e6; % Carrier frequency (Hz)

% Complete these expressions using the variables above
N = Fs*Tmax; % Total number of sample points in the simulation
Ns = Rs*Tmax; % Number of symbols to send
sps = Fs/Rs; % Number of samples per symbol.

% Use these variables for plotting
t = linspace(0, Tmax, N); % Time vector 
f = linspace(-Fs/2, Fs/2, N); % Frequency vector. Recall from 113 that DFT gives us -f/2 to +f/2

%% Symbol generation

% transBits = randi([0, 1], 1, Ns);
% transSymbols = 2*transBits - 1;

bits =[0 1 0 0 0 0 1 1 0 1 1 0 1 1 1 1 0 1 1 0 1 1 1 0 0 1 1 0 0 1 1 1 0 1,...
       1 1 0 0 1 0 0 1 1 0 0 0 0 1 0 1 1 1 0 1 0 0 0 1 1 1 0 0 1 1 0 0 1 0,...
       0 0 0 1 0 0 1 0 0 0 0 0 0 1 0 1 0 0 1 1 0 1 1 0 1 0 0 1 0 1 1 0 1 1,...
       0 1 0 1 1 1 0 1 0 1 0 1 1 0 1 1 0 0 0 1 1 0 0 0 0 1 0 1 1 1 0 1 0 0,...
       0 1 1 0 1 0 0 1 0 1 1 0 1 1 1 1 0 1 1 0 1 1 1 0 0 0 1 0 0 0 0 0 0 1,...
       1 0 1 0 0 1 0 1 1 1 0 0 1 1 0 0 1 0 0 0 0 0 0 1 1 0 0 0 1 1 0 1 1 0,...
       1 1 1 1 0 1 1 0 1 1 0 1 0 1 1 1 0 0 0 0 0 1 1 0 1 1 0 0 0 1 1 0 0 1,...
       0 1 0 1 1 1 0 1 0 0 0 1 1 0 0 1 0 1];

loadStruct = load('module4_symbols.mat');
transSymbols = loadStruct.module4_symbols;

upsampTransSymbols = upsample(transSymbols, sps);

%% Pulse shaping
rrc = rcosdesign(0.2, 10, sps, "sqrt"); 
Xbb = conv(upsampTransSymbols, rrc, "same");

%% Modulation
frequencyOffset = 1e3
complexCarrier = exp(-1j*2*pi*(fc+frequencyOffset)*t);
phaseOffset = pi/4
modulatedSignal = Xbb.*complexCarrier*exp(-1i*phaseOffset);

%% Channel simulation
Tx = real(modulatedSignal);

B = 350e3; %Bandwidth of our WRAP system (f)
s = tf('s'); % Generate the channel filter transfer function
Hc = B*s/(s^2 + s*B + fc^2);
Hd = c2d(Hc, 1/Fs, 'tustin');
[a, b] = tfdata(Hd);

std_ = 0.05
noise = std_ * randn(1,N);
Rx = Tx + noise;
% Filter your transmitted signal with the generated discrete TF. 
% Input variable is "Rx" and output variable is "y"

y = real(filter(a{:}, b{:}, Rx)); %Rx

%% Demodulation

% y = y * 2; % multiply by 2...?

f = [0, 0.1, 0.25, 1];
a = [1, 1, 0, 0];
order = 5; % num coefficients = order + 1
lp = firpm(order, f, a);

ph = zeros(1, N+1);
inph = zeros(1, N);
quad = zeros(1, N);
inph_ = zeros(1, N);
quad_ = zeros(1, N);
err = zeros(1, N);

kp=8.5;
ki=0.1;
integrator=0;

figure;
tl = tiledlayout(2,2,'TileSpacing','Compact');
title(tl, 'Costas loop');
txt = ['PhaseOffset = ' num2str(phaseOffset) ', FreqOffset = ' num2str(frequencyOffset)];
subtitle(tl, txt);
xlabel(tl, 'Index of sample');
ylabel(tl, 'Value of sample');

% y = (y-mean(y))/std(y)/10;%25;
ax1 = nexttile;
y = (y-mean(y))/std(y)/10;
plot(y);
title('Normalized signal');

% inph_(1:order) = y(1:order).*2.*cos(2*pi*fc.*t(1:order)+ph(1:order));
% quad_(1:order) = y(1:order).*-2.*sin(2*pi*fc.*t(1:order)+ph(1:order));

for i=order+1:N
    c = 2*cos(2*pi*fc*t(i)+ph(i));
    s = -2*sin(2*pi*fc*t(i)+ph(i));
    
    inph_(i) = y(i)*c;
    quad_(i) = y(i)*s;

    I = conv(inph_(i-order:i), lp, 'valid');
    Q = conv(quad_(i-order:i), lp, 'valid');
    inph(i) = I(end); % not necessary lol
    quad(i) = Q(end);

    err(i)=inph(i)*quad(i);
    integrator = integrator + ki*err(i);
    ph(i+1)=ph(i)+err(i)*kp+integrator;
end

% inph_ = y.*cos(2*pi*fc*t); 
% quad_ = -1*y.*sin(2*pi*fc*t);
% 
% inph = 2*conv(inph_, lp, 'same');%`filter` doesn't account for the group delay of the filter -- smearing horizontally. what accounts for quadrature??
% quad = 2*conv(quad_, lp, 'same');

%convolving the whole thing -- no delay. but when you do bit by bit,
%there's a delay??



%{
Line up the data. 

tt = tn(1:end-delay); Remove the last `delay` samples of the original and of the time vector.
sn = xn(1:end-delay);

sf = xf;
sf(1:delay) = []; Shift the filtered signal by removing its first `delay` samples.
%}
% 
% t = 0:1/Fs:Tmax-1/Fs; %same as first t vector
% use costas loop to demodulate
% works from ~0.95e6-1.05e6



% df = 0
% f0 = fc+df; %1.05e6; % estimated frequency
% ph = zeros(1, N+1);
% inph = zeros(1, N);
% quad = zeros(1, N);
% inph_ = zeros(1, N);
% quad_ = zeros(1, N);
% err = zeros(1, N);
% % nyquist rate: Fs/2=2000000
% % create lowpass as percent of this
% f = [0, 0.1, 0.15, 1];
% a = [1, 1, 0, 0]; % 0 0
% order = 5;
% lp = firpm(order, f, a);
% integrator = 0;
% 
% kp = 0;
% ki = 0;
% 
% %order=0; %%%%% taking out the filter
% for i = order+1:N
%     c = 2*cos(2*pi*f0*t(i)+ph(i));
%     s = -2*sin(2*pi*f0*t(i)+ph(i));
%     % demodulate the signal by mixing with NCO signals
%     inph_(i) = y(i) * c;
%     quad_(i) = y(i) * s;
%     % low pass filter to get baseband I and Q
%     %shifting back
%     I = filter(lp, 1, inph_(i-order:i));
%     Q = filter(lp, 1, quad_(i-order:i));
%     inph(i) = I(end);%inph_(i);
%     quad(i) = Q(end);%quad_(i);
% 
%     % mix in_ph and quad to get error signal
%     err(i) = inph(i)*quad(i);
% 
%     % loop filter
%     integrator = integrator + ki * err(i);
%     ph(i+1) = ph(i)+kp*err(i)+integrator;
% end




ax2 = nexttile;
hold on;
plot(inph);
plot(quad);
hold off;
title('Inph and quad');
legend('Inph', 'Quad');
% scatter(inph,quad);
% title('Received constellation before SRRC');


% scatter(inph, quad)
% inph(end-5:end)=zeros(1, 6); %shift back
% quad(end-5:end)=zeros(1, 6);
ax3 = nexttile;
plot(1:N, err);
title("Error Plot");

ax4 = nexttile;
plot(1:N+1, ph);
title("Phase Tracker")

linkaxes([ax1 ax2 ax3 ax4], 'x');
xlim([0, 1000]);

reconstructedSignal = inph + 1j*quad;

%% Receive filtering
Ybb = conv(reconstructedSignal, rrc, "same");

sampYbb = Ybb(1:sps:N);
scatterplot(sampYbb);
title('Received constellation');

%% Symbol detection
recSymbols = sign(real(sampYbb));

upsampRecSymbols = upsample(recSymbols, sps);
recBits = 0.5*(recSymbols + 1);

ber = sum(transSymbols ~= recSymbols)/Ns

cum_err = zeros(1, Ns);
cum_err(1) = transSymbols(1) ~= recSymbols(1);
for i=2:Ns
    cum_err(i) = cum_err(i-1)+(transSymbols(i) ~= recSymbols(i));
end

figure;
plot(cum_err);
title("Cumulative bit error");

%% Plots
% Part 2 Transmitted vs Received Signals
figure();
tl = tiledlayout('vertical','TileSpacing','tight');
title(tl, 'Transmitted vs Received Signals');
txt = ['BER = ' num2str(ber) ', Noise = ' num2str(std_)];
subtitle(tl, txt);
xlabel(tl, 'Time (s)');
ylabel(tl, 'Signals');

ax1 = nexttile;
hold on;
plot(t, Xbb*10);
isNZ=(~upsampTransSymbols==0); % addressing logical array of nonzero elements
stem(t(isNZ), upsampTransSymbols(isNZ));
hold off;
title('Transmitted (after first SRRC)');
legend({'Transmitted Baseband * 10', 'Transmitted Symbols'});

ax2 = nexttile;
hold on;
plot(t, Ybb);
isNZ=(~upsampRecSymbols==0); % addressing logical array of nonzero elements
stem(t(isNZ), upsampRecSymbols(isNZ));
hold off;
title('Received (after second SRRC)');
legend({'Received Baseband', 'Received Symbols'});

linkaxes([ax1 ax2],'x');
xlim([0, .001]);