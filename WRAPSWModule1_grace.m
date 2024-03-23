clear all
close all
loadStruct = load('module4_symbols.mat');
transSymbols = loadStruct.module4_symbols;

%part 1
% Values that we're giving you
Fs = 4e6; % Sampling rate (samples/sec)
Rs = 50e3; % Symbol rate in symbols/sec (baud)
Tmax = 0.1; % Max time for the simulation (sec)
fc = 1e6; % Carrier frequency (Hz)

% Complete these expressions using the variables above
N = 4e5; % Total number of sample points in the simulation
Ns = 5000; % Number of symbols to send
sps = 80; % Number of samples per symbol.

% Use these variables for plotting
t = linspace(0 , Tmax,N); % Time vector
f = linspace( -2e6,2e6 , N); % Frequency vector. Recall from 113 that DFT gives us -f/2 to +f/2

B = 350e3; %Bandwidth of our WRAP system (f)
s = tf('s'); % Generate the channel filter transfer function
Hc = B*s/(s^2 + s*B + fc^2);
Hd = c2d(Hc, 1/Fs, 'tustin');
[a, b] = tfdata(Hd);

%binVector = randi([0,1],Ns,1);
bpskVector= transpose(loadStruct.module4_symbols);
upsampledVector = upsample(bpskVector,sps);
h = ones(sps,1);
baseband = conv(upsampledVector,h,"same");
%plot(t,upsampledVector);
%plot(t,baseband);

noise_std = 0.1;
noise = noise_std.*randn(Ns*sps,1);
phaseoffset = 0;
freqoffset = 0;

cc = exp(2 *pi *fc *t *-1i);
%cc = exp(2 *pi *fc *t *-1i-1i*phaseoffset);
%cc = exp(2 *pi *(fc+freqoffset) *t *-1i-1i*phaseoffset);
analytic = baseband.*transpose(cc);

%transmit = real(filter(a{:}, b{:}, analytic+noise));
%transmit = real(analytic)+noise;
transmit = real(analytic);

y = transpose(transmit);
normalization_factor = 2; % adjust this number as needed to recover a nice constellation with points at about -1 and +1
y = (y-mean(y))/std(y)/normalization_factor;

f = [0, 0.1, 0.25, 1];
a = [1, 1, 0, 0];
order = 5;
lp = firpm(order,f,a);
ph = zeros(1, N+1); % "guess" phase for the demodulating sin() and cos()
inph_ = zeros(1, N); % transmitted signal * cos()
quad_ = zeros(1, N); % transmitted signal * sin()
inph = zeros(1, N); % low-pass filtered version of inph_
quad = zeros(1, N); % low-pass filtered version of quad_
err = zeros(1, N); % phase error
inph_(1:order) = y(1:order).*2.*cos(2*pi*fc.*t(1:order)+ph(1:order));
quad_(1:order) = y(1:order).*-2.*sin(2*pi*fc.*t(1:order)+ph(1:order));
integrator = 0;
ki = 0.01;
kp = 0.5;


errsum = 0;
errsps = 0;
phase = -0.0001;
sps_guess = sps;
sps_guess_over_time = zeros(1,N);
for i = 6:N
    inph_(i) = 2*y(i)*cos(2*pi*(fc-115000)*t(i)+ph(i)); 
    quad_(i) = -2*y(i)*sin(2*pi*(fc-115000)*t(i)+ph(i));
    inph(i) = conv(inph_(i-order:i), lp, 'valid');
    quad(i) = conv(quad_(i-order:i), lp, 'valid');
    err(i) = inph(i)*quad(i);
    integrator = integrator + ki*err(i); 
    ph(i+1) = ph(i) + kp*err(i) + integrator; 
end

reconstruct = inph + 1j* quad;
reconstructSamples = zeros(1,N);
% for i = 1:N
%     if(real(reconstruct(i))>0)
%         reconstructSamples(i) = 1;
%     else 
%         reconstructSamples(i) = -1;
%     end
% end 
received = zeros(1,Ns); % we don't know how many symbols were sent -- tx is continuously sending symbols
j = 1;
zeroCrossings = zeros(1,N);
zeroCrossings(1) = 1;
margin = 0.05;
kis = 0.001;
kps = 0.1;
lastZero = 1;
for i = 2:N 
    zeroCrossings(i) = ~(reconstructSamples(i)*reconstructSamples(i-1)+1);
end
sps_guess_over_time(1) = sps_guess;
err_over_time = zeros(1,N);

sps_guess = 70; 
for i = 2:N 
	prev_phase = phase;
	phase = phase + 2*pi/sps_guess;
    %phase = wrapToPi(phase);
    sps_guess_over_time(i) = sps_guess;
    if ( wrapToPi(phase) < -pi*margin && wrapToPi(prev_phase)>pi*margin && i < N-sps_guess)
	   received(j) = reconstructSamples(i);
       j = j+1;
    end

    if (zeroCrossings(i))
        
        distance = i - lastZero;
        lastZero = i;
        modDistance = mod(distance,sps_guess);
        if(modDistance>sps_guess/2)
            errsps = modDistance-sps_guess;
        else
            errsps = modDistance;
        end
        
        %errsps = wrapToPi(phase);
        err_over_time(i) = errsps;
        errsum = errsum + errsps;
	    sps = sps_guess + kps*errsps + kis*errsum;
        sps_guess = sps;
    end
end
plot(err_over_time);
title('error over time');
%figure;
plot(sps_guess_over_time);
sps = round(mean(sps_guess_over_time));
%plot(err)
scatterplot(received);

%k=1;
%%for it = 1:sps_guess:N
 %   received(k) = reconstruct(it);
 %   k = k+1;
%end
%scatterplot(received)
%{
receivedSymbol = zeros(Ns,1);
k = 1;
for i = 1:sps:N
    if real(reconstruct(i)) > 0
        receivedSymbol(k)=1;
    else
        receivedSymbol(k)=-1;
    end
    k = k + 1;
end
scatterplot(receivedSymbol)
%calculating bit error
count = 0;
for i = 1:Ns
    if receivedSymbol(i) ~= bpskVector(i)
        count = count + 1;
    end
end
count;
%test = isequal(receivedSymbol,binVector);
ber = count/Ns;
%noise_std=0;

Xbb = baseband;
Ybb = reconstruct;
upsampTransSymbols = upsampledVector;
upsampRecSymbols = upsample(receivedSymbol,sps);


figure;
tl = tiledlayout(2,2,'TileSpacing','Compact');
title(tl, 'Costas loop validation');
txt = ['PhaseOffset = ' int2str(phaseoffset) ', FreqOffset = ' int2str(freqoffset)];
subtitle(tl, txt);
xlabel(tl, 'Index of sample');
ylabel(tl, 'Value of sample');

ax1 = nexttile;
plot(y./abs(y));
title('Normalized signal'); % don't worry about this for now; basically, it's a sketchy way of doing something called "channel equalization" which we can explore in SW R&D
subtitle('Should be normalized to +/-1');

% <YOUR COSTAS LOOP HERE>

ax2 = nexttile;
hold on;
plot(inph);
plot(quad);
hold off;
title('Inph and quad');
subtitle('Over time, inph should carry most of the signal energy while quad should approach 0'); 
legend('Inph', 'Quad');

ax3 = nexttile;
plot(1:N, err);
title('Error Plot');
subtitle('Should oscillate and then go to zero once the loop locks'); % If your error plot looks like it's at a high frequency (similar to the passband 'Normalized Signal') then you should decrease the upper limit of your LPF. We don't want a control loop whose inputs are oscillating that wildly.

ax4 = nexttile;
plot(1:N+1, ph);
title('Phase Guess');
subtitle('Should initially oscillate and then stabilize once the loop locks');

linkaxes([ax1 ax2 ax3 ax4], 'x');
xlim([0, 5000]); 
%}
