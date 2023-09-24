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

rec_data_bits = run_rx_sim(y, rx_apriori, constants);
disp(rec_data_bits);

% % increment in time
% t = 0:1/Fs:Tmax-1/Fs;
% % use costas loop to demodulate
% % works from ~9e5 to 1.04e6. Locks onto second packet from 8.5e5-1.5
% f0 = 1.04e6; % estimated frequency
% ph = zeros(1, N+1);
% inph = zeros(1, N);
% quad = zeros(1, N);
% inph_ = zeros(1, N);
% quad_ = zeros(1, N);
% err = zeros(1, N);
% % nyquist rate: Fs/2=2500000
% % create lowpass as percent of this
% f = [0, 0.2, 0.25, 1];
% a = [1, 1, 0, 0];
% order = 5;
% lp = firpm(order, f, a);
% integrator = 0;
% 
% kp = 10;
% ki = 0.01;
% 
% for i = order+1:N
%     c = 2*cos(2*pi*f0*t(i)+ph(i));
%     s = -2*sin(2*pi*f0*t(i)+ph(i));
%     % demodulate the signal by mixing with NCO signals
%     inph_(i) = y(i) * c;
%     quad_(i) = y(i) * s;
%     % low pass filter to get rid of high frequency components
%     inph(i-order:i) = filter(lp, 1, inph_(i-order:i));
%     quad(i-order:i)  = filter(lp, 1, quad_(i-order:i));
%     
% 
%     % mix in_ph and quad to get error signal
%     err(i) = inph(i)*quad(i);
% 
%     % NCO
%     integrator = integrator + ki*err(i);
%     ph(i+1) = ph(i)+kp*err(i) + integrator;
% end
% inph(end-5:end)=zeros(1, 6);
% quad(end-5:end) =zeros(1, 6);
% figure();
% plot(1:N, err)
% title("Error Plot");
% 
% figure();
% plot(1:N+1, ph)
% title("Phase Tracker")
% 
% 
% % SRRC
% 
% span = 5;
% sps = 100;
% srrc = rcosdesign(1, span, sps, 'sqrt');
% l = inph+j*quad;
% l = upfirdn(l, srrc,1,1);
% l = l(length(srrc)/2+0.5:end-length(srrc)/2+0.5);
% 
% 
% inph = real(l);
% quad = imag(l);
% figure();
% plot(inph);
% hold on
% plot(quad)
% legend("In_phase", "quad")
% 
% figure();
% scatter(inph, quad)
% title("Constellation Diagram")
% 
% %% Clock recovery
% i = real(l);
% q = imag(l);
% %create zero crossing signal
% zero_cross = ~(sign(i(1:end-1).*i(2:end))+1);
% 
% %Get data rate
% sps = 100;
% 
% %PLL clock recoverypl
% Kp_PLL = 0.3;
% Ki_PLL = 1.1;
% freq = sps;
% phase_acc = pi+0.00001;
% prev_acc = pi-0.00001;
% real_samp = [];
% quad_samp=[];
% err = [];
% integrator = freq;
% for k = 1:N-1
%   %Check if we have crossed zero on the phase accumulator
%   %If the previous value was around -pi and the current value is around pi - we crossed it.
%   if (wrapToPi(phase_acc )< -pi * 0.75 && wrapToPi(prev_acc) > pi*0.75)
%     %we have crossed zero - sample
%     real_samp(end+1) = i(k);
%     quad_samp(end+1) = q(k);
%   end
%   %Detect zero crossing and use it to tune PLL
%   if (zero_cross(k))
%     %we have a zero crossing - adjust freq
%     err(end+1) = wrapToPi(phase_acc);
%     integrator = integrator + err(end)*Ki_PLL;
%     freq = err(end) * Kp_PLL + integrator;
%   end
%   prev_acc = phase_acc;
%   phase_acc = phase_acc + 2*pi/freq;
% end
% 
% stem (real_samp);
% plot(quad_samp);
% selected = real_samp +j*quad_samp;
% legend("Real sampled", "Imaginary sampled");
% figure;
% scatterplot(real_samp+quad_samp*j);
% figure;
% plot(err);
% legend("PLL error");
% 
% %% Bit conversion
% %Conversion to bits
% selected_sign = abs(real(selected))./real(selected);
% bits_received = double(selected_sign > 0);
% key2 = [1, 1, 1, -1, -1, -1, 1, -1, -1, 1, 1, -1, 1, -1,1];
% [bit_corr, bit_pos] = xcorr( (bits_received -0.5)*2, key2);
% plot(bit_pos, bit_corr);
% legend("Bit correlation to key");
% shift = 0;
% 
% for i=1:size(bit_corr,2)
%     if abs(bit_corr(i) + bit_corr(i+15) +bit_corr(i+30) ) > 39
% 
%         shift = bit_pos(i)+45;
%         if (bit_corr(i) < 0)
%             bits_received = bits_received * -1 + 1;
%         end
%         break;
%     end
% end
% received = bits_received(shift+1:shift+256);
% strarr = int2str(received);
% disp(strarr);
%}


