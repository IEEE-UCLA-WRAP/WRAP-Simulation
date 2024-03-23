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

% bits =[0 1 0 0 0 0 1 1 0 1 1 0 1 1 1 1 0 1 1 0 1 1 1 0 0 1 1 0 0 1 1 1 0 1,...
%        1 1 0 0 1 0 0 1 1 0 0 0 0 1 0 1 1 1 0 1 0 0 0 1 1 1 0 0 1 1 0 0 1 0,...
%        0 0 0 1 0 0 1 0 0 0 0 0 0 1 0 1 0 0 1 1 0 1 1 0 1 0 0 1 0 1 1 0 1 1,...
%        0 1 0 1 1 1 0 1 0 1 0 1 1 0 1 1 0 0 0 1 1 0 0 0 0 1 0 1 1 1 0 1 0 0,...
%        0 1 1 0 1 0 0 1 0 1 1 0 1 1 1 1 0 1 1 0 1 1 1 0 0 0 1 0 0 0 0 0 0 1,...
%        1 0 1 0 0 1 0 1 1 1 0 0 1 1 0 0 1 0 0 0 0 0 0 1 1 0 0 0 1 1 0 1 1 0,...
%        1 1 1 1 0 1 1 0 1 1 0 1 0 1 1 1 0 0 0 0 0 1 1 0 1 1 0 0 0 1 1 0 0 1,...
%        0 1 0 1 1 1 0 1 0 0 0 1 1 0 0 1 0 1];

loadStruct = load('module4_symbols.mat');
transSymbols = loadStruct.module4_symbols;

% num_sample_delay = 1;%round(normrnd(500, 100));
% sample_delay = zeros(1, num_sample_delay);
% transSymbols = [sample_delay transSymbols(1:Ns-num_sample_delay)];

upsampTransSymbols = upsample(transSymbols, sps);

%% Pulse shaping
rrc = rcosdesign(0.2, 10, sps, "sqrt"); 
Xbb = conv(upsampTransSymbols, rrc, "same");

%% Delay for symbol sync
% FD = 0.32381; % fractional delay
% num_taps = 100; % # of filter taps (coefficients) 
% % h = designFracDelayFIR(FD,num_taps);
% % fdfir = dsp.FIRFilter(h); % magic Matlab function
% test = sin(t);
% % y = fdfir(test);
% d = fdesign.fracdelay((1/Fs)*0.5,Fs);
% fdfir = design(d, 'lagrange', 'FilterStructure', 'farrowfd');
% y = fdfir(test);
% plot(t,y, t,test);
% legend('Filter Output','Original Sequence');
% title('Raw Filter Output v.s. Input Sequence');

% delay = 30; % fractional delay, in samples 60.4->12, 12
% % N_delay = 21; % number of taps
% % n = round(-N_delay/2):round(N_delay/2); % ...-3,-2,-1,0,1,2,3...
% % h = sinc(n - delay); % # calc filter taps
% % h = h .* hamming(N_delay); %# window the filter to make sure it decays to 0 on both sides
% % h = h/sum(h); % # normalize to get unity gain, we don't want to change the amplitude/power
% delayed_Xbb = [zeros(1,delay) Xbb(1:N-delay)];%conv(Xbb, h, 'full'); %# apply filter
% figure;
% hold on;
% stem(Xbb(1000:1100));
% stem(delayed_Xbb(1000:1100));
% % delayed_Xbb = delayed_Xbb(1:N);
% %Xbb = delayed_Xbb;
% hold off;
% legend('Original Sequence', 'Filter Output');
% title('Raw Filter Output v.s. Input Sequence');

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

%%
%%%%%%%%%%%%%%%%%%%%%%%%%% this is different than frame sync because it
%%%%%%%%%%%%%%%%%%%%%%%%%% doesn't solve false lock
% Next, let’s add a random channel propagation delay in units of sampling intervals (not symbol intervals):
timeOffset = 0; % Delay (in samples) added [20->late by 20, 160->early by 20) 20->err
% Delayed sequence
y = [zeros(1, timeOffset), Rx(1:end-timeOffset)];
%%%%%%%%%%%%%%%%%%%%%%%%%%


% y = Rx;%real(filter(a{:}, b{:}, Rx)); %Rx

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
CL_errs = zeros(1, N);

kp=8.5;
ki=0.1;
integrator=0;

figure;
tl = tiledlayout(2,2,'TileSpacing','Tight');
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
    inph(i) = I;
    quad(i) = Q;

    CL_errs(i)=inph(i)*quad(i);
    integrator = integrator + ki*CL_errs(i);
    ph(i+1)=ph(i)+CL_errs(i)*kp+integrator;
end

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
plot(1:N, CL_errs);
title("Error Plot");

ax4 = nexttile;
plot(1:N+1, ph);
title("Phase Tracker")

linkaxes([ax1 ax2 ax3 ax4], 'x');
xlim([0, 1000]);

reconstructedSignal = inph + 1j*quad;

%% Receive filtering
Ybb = conv(reconstructedSignal, rrc, "same");
save("symb_sync_input.mat","Ybb");
% sps;% Number of samples per symbol.

%% Receiver downsampling

% sps_receive = sps;%(Fs+1)/Rs; % instead of Fs/Rs, the "true" sps of the transmitter
% rx_clock_offset = 1; % # of samples in "offset" between the tx clock and rx clock
% sampYbb = Ybb(rx_clock_offset:sps_receive:N-rx_clock_offset); % instead of Ybb(1:sps:N)
% sampYbb = Ybb(1:sps:N);

inph_Ybb = real(Ybb);
quad_Ybb = imag(Ybb);
%create zero crossing signal
zero_cross = ~(sign(inph_Ybb(1:end-1).*inph_Ybb(2:end))+1);

%% Alternative approach
% received = zeros(1,Ns);
% j = 1;
% zeroCrossings = zeros(1,N);
% zeroCrossings(1) = 1;
% margin = 0.05;
% kis = 0;%0.001;
% kps = 0;%0.1;
% lastZero = 1;
% % for i = 2:N 
% %     zeroCrossings(i) = ~(reconstructSamples(i)*reconstructSamples(i-1)+1);
% % end
% zeroCrossings = ~(sign(inph_Ybb(1:end-1).*inph_Ybb(2:end))+1);
% sps_guess = 80;
% sps_guess_over_time(1) = sps_guess;
% err_over_time = zeros(1,N);
% for i = 2:N-1 
% 	% prev_phase = phase;
% 	% phase = phase + 2*pi/sps_guess;
%     % %phase = wrapToPi(phase);
%     % sps_guess_over_time(i) = sps_guess;
%     % if ( wrapToPi(phase) < -pi*margin && wrapToPi(prev_phase)>pi*margin && i < N-sps_guess)
% 	%    received(j) = inph_Ybb(i);
%     %    j = j+1;
%     % end
% 
%     if (zeroCrossings(i))
% 
%         distance = i - lastZero;
%         lastZero = i;
%         modDistance = mod(distance,sps_guess);
%         if(modDistance>sps_guess/2)
%             errsps = modDistance-sps_guess;
%         else
%             errsps = modDistance;
%         end
% 
%         %errsps = wrapToPi(phase);
%         err_over_time(i) = errsps;
%         % errsum = errsum + errsps;
% 	    % sps = sps_guess + kps*errsps + kis*errsum;
%         % sps_guess = sps;
%     end
% end
% 
% figure; plot(err_over_time);
%% 

window = 100;
coeffs = ones(1, window)/window;

sps_PLL=81;
Kp_PLL = 0;%1.0;%0.3;
Ki_PLL = 0;%0.01;%1.1;

real_samp = [];
quad_samp=[];
samp_index = [];
ZC_errs = [];
% ZC_mod_errs = [];
smoothed_errs = [];
smoothed_errs_idx = [];

sps_guess = zeros(1,N-1);
sps_guess(1) = sps_PLL;

phase_acc = pi+0.00001;
prev_acc = pi-0.00001;
phase_PLL = zeros([1, N-1]);
phase_PLL(1) = prev_acc;
phase_acc_array = [];

integrator = 0;
err_window = 0;
% lastZero = 1;

for k = 2:N-1
  sps_guess(k) = sps_guess(k-1); % prev unless it gets changed in the if-if
  %Check if we have crossed zero on the phase accumulator
  %If the previous value was around -pi and the current value is around pi - we crossed it.
  if (wrapToPi(phase_acc )< -pi * 0.75 && wrapToPi(prev_acc) > pi*0.75)
    %we have crossed zero - sample
    real_samp(end+1) = inph_Ybb(k); % will change every iteration but dw about it -- not sure if we do this in the C code
    quad_samp(end+1) = quad_Ybb(k);
    samp_index(end+1) = k;
  end
  %Detect zero crossing and use it to tune PLL
  if (zero_cross(k))
    %we have a zero crossing - adjust freq
    err_window = err_window + 1;
    ZC_errs(end+1) = wrapToPi(phase_acc); % not directly a function of sps_guess

    if (length(ZC_errs) == 99)
       disp("stop")
    end

    % if (mod(err_window,window) == 0) %including this introduces NaN's?!?!
    %     smoothed_err = sum(ZC_errs(end-99:end))/(window+1);
    %     smoothed_errs(end+1) = smoothed_err;
    %     smoothed_errs_idx(end+1) = k;
    %     integrator = integrator + smoothed_err*Ki_PLL; % err(end)
    %     sps_PLL = smoothed_err * Kp_PLL + integrator; % err(end)
    %     sps_guess(k) = sps_guess(k) + sps_PLL;
    % end
    
  end

  prev_acc = phase_acc;
  phase_acc = phase_acc + 2*pi/sps_PLL;
  % phase_acc_array(end+1) = phase_acc;

  phase_PLL(k) = wrapToPi(phase_acc);

end

length(samp_index)

figure();
plot (real_samp);
hold on;
plot(quad_samp);
hold off;
sampYbb_PLL = real_samp +1j*quad_samp;
legend("Real sampled", "Imaginary sampled");

scatterplot(real_samp+quad_samp*1j);
title('Received constellation');

% (BEGIN) Symbol sync error function validation
figure;
tl = tiledlayout(2,1, 'TileSpacing', 'tight');
title(tl, 'Symbol sync error function validation');
subtitle(tl, 'Zero-crossings do not occur EXACTLY in between symbols');
xlabel(tl, 'Sample #');
txt = ['Transmitter sps: ' num2str(sps) ', K_p PLL = ' num2str(Kp_PLL) ', K_i PLL = ' num2str(Ki_PLL)];
subtitle(tl, txt);

ax1 = nexttile;
hold on;
p1 = stem(1:sps:N, upsampTransSymbols(1:sps:N),'Color',"#0072BD"); %default blue
ideal_zero_cross = ~(sign(transSymbols(1:end-1).*transSymbols(2:end))+1);
ideal_zero_cross = [zeros(1,40) upsample(ideal_zero_cross, sps)];
ideal_zc_indeces = find(ideal_zero_cross);

for k=1:10
   p2 = xline(ideal_zc_indeces(k),'Color',"#7E2F8E",'Label',ideal_zc_indeces(k),'LabelOrientation','horizontal','LineWidth',3);
end
p3 = stem(inph_Ybb, 'LineStyle','none', 'Color',"#EDB120"); % default orange

hold off;
legend([p1 p2 p3],'Upsampled TRANSMITTED symbols (perfectly upsampled by 80)', 'TX-INTENDED zero-crossing locations (at perfect sps/2 intervals)', 'RECEIVED baseband signal');
% title('Zero-crossings do not occur EXACTLY in between symbols');
grid on;
grid minor;

ax2 = nexttile;
hold on;
for k=1:10
    xline(ideal_zc_indeces(k),'Color',"#7E2F8E",'Label',ideal_zc_indeces(k),'LabelOrientation','horizontal','LineWidth',3);
end
p1 = stem(phase_PLL, 'LineStyle','none', 'Color', "#EDB120"); % default orange
zc_indeces = find(zero_cross); % reusing var name
p2 = stem(zc_indeces,ZC_errs, 'LineWidth',3, 'Color',"#7E2F8E"); % default purple
% p3 = stem(smoothed_errs_idx, smoothed_errs,'LineWidth',3);
hold off;
% title('Symbol sync validation');
legend([p1 p2],'PLL phase (its frequency is controlled by sps-PLL)', 'PLL error (at locations of the RX-DETECTED zero-crossings)');
for k=1:10
      text(zc_indeces(k)+5,ZC_errs(k),['(' num2str(zc_indeces(k)) ',' num2str(ZC_errs(k)) ')'],'Color',"#7E2F8E");
end

linkaxes([ax1 ax2], 'x');
xlim([1 800]);
grid on;
grid minor;

% (END) Symbol sync error function validation

figure;
tl = tiledlayout(2,1, 'TileSpacing', 'tight');
title(tl, 'Timing Recovery K_P,K_I tuning')
ax1 = nexttile;
hold on;
stem(find(zero_cross), ZC_errs); % should be close to 0, plot every time there's a zero-crossing
stem(smoothed_errs_idx, smoothed_errs, 'LineWidth',3);
hold off;
title(ax1,'PLL error over time');
legend(ax1,"PLL error", "PLL smoothed error");
xlabel(ax1,'Sample #');
grid on;
grid minor;

ax2=nexttile;
stem(sps_guess, 'LineStyle','none');
title(ax2,'SPS guess over time');
xlabel(ax2,'Sample #')

linkaxes([ax1 ax2], 'x');
% xlim([1 1e5]);
grid on;
grid minor;

%%

% scatterplot(sampYbb);

%Symbol detection
recSymbols = sign(real(sampYbb_PLL));

% upsampRecSymbols = upsample(recSymbols, sps);
% upsampRecSamples = upsample(sampYbb, sps_receive);
recBits = 0.5*(recSymbols + 1);

min_symbs = min(length(transSymbols),length(recSymbols))
ber = sum(transSymbols(1:min_symbs) ~= recSymbols(1:min_symbs))/Ns

cum_err = zeros(1, min_symbs);
cum_err(1) = transSymbols(1) ~= recSymbols(1);
for i=2:min_symbs
    cum_err(i) = cum_err(i-1)+(transSymbols(i) ~= recSymbols(i));
end
first_error_index = find(cum_err,1);

figure;
plot(cum_err, 'LineWidth',3);
title("Cumulative bit error");
xlabel('Sample #');
ylabel('# of cumulative bit errors');
subtitle([" sps @ receiver: " + num2str(sps_PLL) " sample index of first error: " + num2str(first_error_index)]);
% xlim([0 100]);
% ylim([-1 30]);

%% Plots
% Part 2 Transmitted vs Received Signals
figure();
tl = tiledlayout('vertical','TileSpacing','tight');
title(tl, 'Transmitted vs Received Signals');
txt = ['BER = ' num2str(ber) ', Noise = ' num2str(std_)];
subtitle(tl, txt);
xlabel(tl, 'Sample (#)');
ylabel(tl, 'Signals');

ax1 = nexttile;
hold on;
plot(Xbb, 'LineWidth',2); %
% isNZ=(~upsampTransSymbols==0); % addressing logical array of nonzero elements
stem((1:sps:N),real(Xbb(1:sps:N)),'LineWidth',2); %t(isNZ), upsampTransSymbols(isNZ)
hold off;
title('Transmitted (after first SRRC)');
legend({'Transmitted Baseband (not *10)', 'Transmitted SAMPLES'});
grid minor;

ax2 = nexttile;
hold on;
plot(real(Ybb),'LineWidth',2);
% isNZ=(~upsampRecSymbols==0); % addressing logical array of nonzero elements
stem(samp_index,real(sampYbb_PLL), 'LineWidth', 2); %t(isNZ), upsampRecSymbols(isNZ) 1:sps:N
hold off;
title('Received (after second SRRC)');
legend({'Received Baseband', 'Received SAMPLES'});
grid minor;

linkaxes([ax1 ax2],'x');
xlim([15e3, 16e3]);