clear all;
close all;
%==========================================
% Title:  Link Simulation
% Authors: David Zheng & Nathan Pereira
% Last Modified: December 10, 2022
%==========================================
% This version of the code modulates and demodulates the signal.
% The packet header is found. 
%%
Fs = 5000000;             % Sampling rate
N  = 300000;               % Number of points sampled
Tmax = N/Fs;              % length of packet (s)
Rs = 50000;               % Symbol rate in symbols/sec (baud)
Ns = Rs * Tmax;           % Number of symbols to send
sps = Fs / Rs;            % Number of samples per symbol
fc = 1e6;             % Carrier frequency in Hz
B = 350000;         % Bandwidth

% for plotting
t = linspace(0, Tmax - 1/Fs, N); % Time variable
f = linspace(-Fs/2,Fs/2, N);        % Frequency variable

% RRC setup
span = 5;           % number of symbols to truncate RRC
rolloff = 0;      % parameter of RRC
RRC = rcosdesign(rolloff, span, sps, 'sqrt');

%% Symbol Generation
% Generate random set of bits
bits_per_packet = 256;

bits = [0 1 0 0 0 0 1 1 0 1 1 0 1 1 1 1 0 1 1 0 1 1 1 0 0 1 1 0 0 1 1 1 0 1,...
        1 1 0 0 1 0 0 1 1 0 0 0 0 1 0 1 1 1 0 1 0 0 0 1 1 1 0 0 1 1 0 0 1 0,...
        0 0 0 1 0 0 1 0 0 0 0 0 0 1 0 1 0 0 1 1 0 1 1 0 1 0 0 1 0 1 1 0 1 1,...
        0 1 0 1 1 1 0 1 0 1 0 1 1 0 1 1 0 0 0 1 1 0 0 0 0 1 0 1 1 1 0 1 0 0,...
        0 1 1 0 1 0 0 1 0 1 1 0 1 1 1 1 0 1 1 0 1 1 1 0 0 0 1 0 0 0 0 0 0 1,...
        1 0 1 0 0 1 0 1 1 1 0 0 1 1 0 0 1 0 0 0 0 0 0 1 1 0 0 0 1 1 0 1 1 0,...
        1 1 1 1 0 1 1 0 1 1 0 1 0 1 1 1 0 0 0 0 0 1 1 0 1 1 0 0 0 1 1 0 0 1,...
        0 1 0 1 1 1 0 1 0 0 0 1 1 0 0 1 0 1];

%{
bits = [0  1  0  0  1  0  0  1  0  0  1  0  0  0  0  0  0  1  1  0  0  1  0,...
        0  0  1  1  0  1  1  1  1  0  1  1  0  1  1  1  0  0  0  1  0  0  1,...
        1  1  0  1  1  1  0  1  0  0  0  0  1  0  0  0  0  0  0  1  1  0  1,...
        0  1  1  0  1  1  0  1  1  1  0  0  1  1  0  1  1  1  1  0  1  1  1,...
        0  1  1  1  0  0  1  0  0  0  0  0  0  1  1  1  1  0  0  1  0  1  1,...
        0  1  1  1  1  0  1  1  1  0  1  0  1  0  0  1  0  0  0  0  0  0  1,...
        1  1  0  1  0  0  0  1  1  0  0  1  0  1  0  1  1  0  1  1  0  0  0,...
        1  1  0  1  1  0  0  0  0  1  0  0  0  0  0  0  1  1  0  1  1  0  1,...
        0  1  1  0  0  1  0  1  0  0  1  0  0  0  0  0  0  1  1  0  0  1  0,...
        1  0  1  1  0  0  1  1  1  0  1  1  0  0  1  1  1  0  1  1  0  0  1,...
        0  1  0  1  1  1  0  0  1  0  0  1  1  1  0  1  0  0  0  0  1  0  0,...
        0  0  1];
%}
% Map bits to BPSK symbols
symbs = 2 * (bits - 0.5);

% create packet key
key = [1, 1, 1, -1, -1, -1, 1, -1, -1, 1, 1, -1, 1, -1, 1];
packet_header = key;
for i=1:2
    packet_header = cat(2, key, packet_header);
end

% add 3 keys to the front of the data. the packet header. 
symbs = cat(2, packet_header, symbs);
p_len = length(symbs);
samples = zeros(1, N/sps);
for i = 1:N/sps
    samples(i) = symbs(mod(i, p_len)+1);
end

% upsample
deltas = upsample(samples,sps);

%% Pulse Shaping (RRC)
symbol_wave = upfirdn(deltas, RRC, 1, 1);
symbol_wave = symbol_wave(length(RRC)/2+0.5:end-length(RRC)/2+0.5);


%% Carrier Modulation
comp_carr = exp(1i*2*pi*fc*t);
analytic_sig = symbol_wave .* comp_carr;
phase = unifrnd(0, 2*pi);
Tx = real(analytic_sig*exp(1i*phase));

%% Noise and Band Limited Channel
std_ = 0.01;
noise = std_ * randn(1,N);
Rx = Tx + noise;

s = tf('s');
Hc = B*s/(s^2 + s*B + fc^2);
Hd = c2d(Hc, 1/Fs, 'tustin');
[a, b] = tfdata(Hd);
  
y = real(filter(a{:}, b{:}, Rx));

%figure()
%plot(t(10000:11000), y(10000:11000));

y = (y-mean(y))/std(y)/25;
figure();
plot(y);
t = 0:1/Fs:Tmax-1/Fs;
% use costas loop to demodulate
% works from ~0.95e6-1.05e6
df = 0;
f0 = fc+df; %1.05e6; % estimated frequency
ph = zeros(1, N+1);
inph = zeros(1, N);
quad = zeros(1, N);
inph_ = zeros(1, N);
quad_ = zeros(1, N);
err = zeros(1, N);
% nyquist rate: Fs/2=2000000
% create lowpass as percent of this
f = [0, 0.1, 0.15, 1];
a = [1, 1, 0, 0];
order = 5;
lp = firpm(order, f, a);
integrator = 0;

kp = 8.5;
ki = 0.1;

for i = order+1:N
    c = 2*cos(2*pi*f0*t(i)+ph(i));
    s = -2*sin(2*pi*f0*t(i)+ph(i));
    % demodulate the signal by mixing with NCO signals
    inph_(i) = y(i) * c;
    quad_(i) = y(i) * s;
    % low pass filter to get rid of high frequency components
    I = filter(lp, 1, inph_(i-order:i));
    Q = filter(lp, 1, quad_(i-order:i));
    inph(i) = I(end);
    quad(i) = Q(end);

    % mix in_ph and quad to get error signal
    err(i) = inph(i)*quad(i);

    % NCO
    integrator = integrator + ki * err(i);
    ph(i+1) = ph(i)+kp*err(i)+integrator;
end
figure()
plot(inph(1:10000))
figure()
scatter(inph, quad)
inph(end-5:end)=zeros(1, 6);
quad(end-5:end)=zeros(1, 6);
figure();
plot(1:N, err)
title("Error Plot");

figure();
plot(1:N+1, ph)
title("Phase Tracker")
A = [1:N+1; ones(length(ph), 1)'];
x = A'\ph';

% SRRC
l = inph+1j*quad;
l = upfirdn(l, RRC,1,1);
l = l(length(RRC)/2+0.5:end-length(RRC)/2+0.5);

inph = real(l);
quad = imag(l);
figure();
plot(inph);
hold on
plot(quad)
hold off
legend("In_phase", "quad")

figure();
scatter(inph, quad)
title("Constellation Diagram")

%% Clock recovery
i = real(l);
q = imag(l);
%create zero crossing signal
zero_cross = ~(sign(i(1:end-1).*i(2:end))+1);

%PLL clock recovery pll
Kp_PLL = 0.3;
Ki_PLL = 1.1;
phase_acc = pi+0.00001;
prev_acc = pi-0.00001;
real_samp = [];
quad_samp=[];
err = [];
integrator = sps;
for k = 1:N-1
  %Check if we have crossed zero on the phase accumulator
  %If the previous value was around -pi and the current value is around pi - we crossed it.
  if (wrapToPi(phase_acc )< -pi * 0.75 && wrapToPi(prev_acc) > pi*0.75)
    %we have crossed zero - sample
    real_samp(end+1) = i(k);
    quad_samp(end+1) = q(k);
  end
  %Detect zero crossing and use it to tune PLL
  if (zero_cross(k))
    %we have a zero crossing - adjust freq
    err(end+1) = wrapToPi(phase_acc);
    integrator = integrator + err(end)*Ki_PLL;
    sps = err(end) * Kp_PLL + integrator;
  end
  prev_acc = phase_acc;
  phase_acc = phase_acc + 2*pi/sps;
end

figure();
plot (real_samp);
hold on;
plot(quad_samp);
hold off;
selected = real_samp +1j*quad_samp;
legend("Real sampled", "Imaginary sampled");
scatterplot(real_samp+quad_samp*1j);
figure;
plot(err);
legend("PLL error");

%% Bit conversion
%Conversion to bits
selected_sign = abs(real(selected))./real(selected);
bits_received = double(selected_sign > 0);
key2 = [1, 1, 1, -1, -1, -1, 1, -1, -1, 1, 1, -1, 1, -1, 1];
[bit_corr, bit_pos] = xcorr( (bits_received -0.5)*2, key2);
plot(bit_pos, bit_corr);
legend("Bit correlation to key");
shift = 0;

for i=1:size(bit_corr,2)
    if abs(bit_corr(i) + bit_corr(i+15) +bit_corr(i+30) ) > 39
        shift = bit_pos(i)+45;
        if (bit_corr(i) < 0)
            bits_received = bits_received * -1 + 1;
        end
        break;
    end
end
received = bits_received(shift+1:shift+256);
strarr = int2str(received);
disp(strarr);

