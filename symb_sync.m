close all;
clear all;

% Values that we're giving you
Fs = 4e6; % Sampling rate (samples/sec)
Rs = 50e3; % Symbol rate in symbols/sec (baud)
Tmax = 0.1; % Max time for the simulation (sec)
fc = 1e6; % Carrier frequency (Hz)

% Complete these expressions using the variables above
N = Fs*Tmax; % Total number of sample points in the simulation
Ns = Rs*Tmax; % Number of symbols to send
sps = N/Ns; % Number of samples per symbol.

loadStruct = load('symb_sync_input.mat');
Ybb = loadStruct.Ybb;

loadStruct = load('module4_symbols.mat');
transSymbols = loadStruct.module4_symbols;

inph_Ybb = real(Ybb);
quad_Ybb = imag(Ybb);

%%
%create zero crossing signal
zero_cross = ~(sign(inph_Ybb(1:end-1).*inph_Ybb(2:end))+1);

%%
figure;
plot(real(Ybb));
title('Symbol sync input (after Costas loop)');
xlim([0 10000]);

scatterplot(Ybb(1:sps:N));
title('IDEAL Received constellation (sampled perfectly)');

%% Alternative approach
% received = zeros(1,Ns);
% j = 1;
% 
% margin = 0.05;
% kis = 0.001;
% kps = 0.1;
% lastZero = 1;
% 
% phase = -0.0001;
% errsum = 0;
% 
% sps_guess = 70;
% sps_guess_over_time(1) = sps_guess;
% err_over_time = zeros(1,N);
% for i = 2:N-1 
% 	prev_phase = phase;
% 	phase = phase + 2*pi/sps_guess;
%     %phase = wrapToPi(phase);
%     sps_guess_over_time(i) = sps_guess;
%     if ( wrapToPi(phase) < -pi*margin && wrapToPi(prev_phase)>pi*margin && i < N-sps_guess)
% 	   received(j) = inph_Ybb(i);
%        j = j+1;
%     end
% 
%     if (zero_cross(i))
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
%         errsum = errsum + errsps;
% 	    sps = sps_guess + kps*errsps + kis*errsum;
%         sps_guess = sps;
%     end
% end
% 
% figure; plot(err_over_time);
%% 

% window = 50;
% coeffs = ones(1, window)/window;

sps_PLL=80;
Kp_PLL = 0.0;%1.0;%0.3;
Ki_PLL = 0.0;%0.01;%1.1;

sampPhaseOffset = 0;
Ybb = Ybb(1+sampPhaseOffset:end);
% sampled_baseband = Ybb(1+sampPhaseOffset:sps:N);
% sampled_symbols = Ybb(1:sps-1:N); % sps-1 or sps-2
% sampled_symbols = sampled_symbols(1:Ns); % Recall that sps=N/Ns. Since your sps has decreased and your N has remained constant (you're still receiving the same # of samples), your Ns will *increase*. So you will end up receiving more than Ns=5000 symbols! Truncate to the first Ns symbols to compare these to the Ns transmitted symbols, to compute your bit error rate.
% sampled_symbols = Ybb(1:sps+1:N); % sps+1 or sps+2
% extra_random_symbols = randi([-1, 1], [1, Ns-numel(sampled_symbols)]); % Recall that sps=N/Ns. Since your sps has increased and your N has remained constant (you're still receiving the same # of samples), your Ns will *decrease*. So you will end up receiving fewer than Ns=5000 symbols! To handle this, just fill the remaining symbols with random symbols, which simulates just picking up ambient noise after the transmitter's done transmitting.
% sampled_symbols = [sampled_symbols extra_random_symbols];

% figure;
% for i=1:length(sampled_symbols)/3
%     scatter(real(sampled_symbols(1:i)), imag(sampled_symbols(1:i)));
%     xlim([-2 2]);
%     ylim([-2 2]);
%     drawnow
%     pause(0.001);
% end
% title('Animated plot of first Ns/4 samples @ receiver');
% scatterplot(sampled_baseband);
% receivedSymbs = sign(real(sampled_baseband));
% BER = sum(transSymbols(1:end-1) ~= receivedSymbs)/Ns;
% title(['Sampling phase offset =' num2str(sampPhaseOffset) ' BER=' num2str(BER)]);

real_samp = [];
quad_samp=[];
samp_index = [];
ZC_errs = [];

sps_guess = zeros(1,N-1);
sps_guess(1) = sps_PLL;

phase_acc = pi+0.00001;
prev_acc = pi-0.00001;
phase_PLL = zeros([1, N-1]);
phase_PLL(1) = prev_acc;

integrator = 0;
lastZero = find(zero_cross,1);

for k = 2:N-1
 
  %Check if we have crossed zero on the phase accumulator
  if (wrapToPi(phase_acc )< -pi * 0.75 && wrapToPi(prev_acc) > pi*0.75)
    %we have crossed zero - sample
    real_samp(end+1) = inph_Ybb(k);
    quad_samp(end+1) = quad_Ybb(k);
    samp_index(end+1) = k;
  end
  %Detect zero crossing and use it to tune PLL
  if (zero_cross(k))
    %we have a zero crossing - adjust freq
    % ZC_errs(end+1) = wrapToPi(phase_acc); % not directly a function of sps_guess
    distance = k - lastZero;
    lastZero = k;
    modDistance = mod(distance,sps_guess(k));
    if(modDistance>sps_guess(k)/2)
        errsps = modDistance-sps_guess(k);
    else
        errsps = modDistance;
    end
    
    ZC_errs(end+1) = errsps;
    errsum = integrator + errsps;
    sps = sps_guess(k-1) + Kp_PLL*errsps + Ki_PLL*integrator;
    sps_guess(k) = sps;
  else % no zero crossing
      sps_guess(k) = sps_guess(k-1); % prev unless it gets changed in the if-if
  end

  prev_acc = phase_acc;
  phase_acc = phase_acc + 2*pi/sps_guess(k);
  phase_PLL(k) = wrapToPi(phase_acc);

end

figure();

plot(real_samp);
hold on;
plot(quad_samp);
hold off;
sampYbb_PLL = real_samp +1j*quad_samp;
title('Sampled values');
legend("Real sampled", "Imaginary sampled");

% Timing Recovery K_P,K_I tuning
figure;
tl = tiledlayout(2,1, 'TileSpacing', 'tight');
title(tl, 'Timing Recovery K_P,K_I tuning')
ax1 = nexttile;
hold on;
stem(find(zero_cross), ZC_errs); % should be close to 0, plot every time there's a zero-crossing
% stem(smoothed_errs_idx, smoothed_errs, 'LineWidth',3);
hold off;
title(ax1,'PLL error over time');
legend(ax1,"PLL error");
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


scatterplot(real_samp+quad_samp*1j);
title('Received constellation');

%% (BEGIN) Symbol sync error function validation
figure;
% tl = tiledlayout(2,1, 'TileSpacing', 'tight');
title(tl, 'Symbol sync error function validation');
subtitle(tl, 'Zero-crossings do not occur EXACTLY in between symbols');
xlabel(tl, 'Sample #');
txt = ['Transmitter sps: ' num2str(sps) ', K_p PLL = ' num2str(Kp_PLL) ', K_i PLL = ' num2str(Ki_PLL)];
subtitle(tl, txt);
% 
% ax1 = nexttile;
% hold on;
% p1 = stem(1:sps:N, upsampTransSymbols(1:sps:N),'Color',"#0072BD"); %default blue
% ideal_zero_cross = ~(sign(transSymbols(1:end-1).*transSymbols(2:end))+1);
% ideal_zero_cross = [zeros(1,40) upsample(ideal_zero_cross, sps)];
% ideal_zc_indeces = find(ideal_zero_cross);
% 
% for k=1:10
%    p2 = xline(ideal_zc_indeces(k),'Color',"#7E2F8E",'Label',ideal_zc_indeces(k),'LabelOrientation','horizontal','LineWidth',3);
% end
% p3 = stem(inph_Ybb, 'LineStyle','none', 'Color',"#EDB120"); % default orange
% 
% hold off;
% legend([p1 p2 p3],'Upsampled TRANSMITTED symbols (perfectly upsampled by 80)', 'TX-INTENDED zero-crossing locations (at perfect sps/2 intervals)', 'RECEIVED baseband signal');
% % title('Zero-crossings do not occur EXACTLY in between symbols');
% grid on;
% grid minor;
% 
% ax2 = nexttile;
% hold on;
% for k=1:10
%     xline(ideal_zc_indeces(k),'Color',"#7E2F8E",'Label',ideal_zc_indeces(k),'LabelOrientation','horizontal','LineWidth',3);
% end
p1 = stem(phase_PLL, 'LineStyle','none', 'Color', "#EDB120"); % default orange
zc_indeces = find(zero_cross); % reusing var name
p2 = stem(zc_indeces,ZC_errs, 'LineWidth',3, 'Color',"#7E2F8E"); % default purple
% p3 = stem(smoothed_errs_idx, smoothed_errs,'LineWidth',3);
hold off;
% title('Symbol sync validation');
legend('PLL error (at locations of the RX-DETECTED zero-crossings)'); % [p1 p2],'PLL phase (its frequency is controlled by sps-PLL)', 
for k=1:10
      text(zc_indeces(k)+5,ZC_errs(k),['(' num2str(zc_indeces(k)) ',' num2str(ZC_errs(k)) ')'],'Color',"#7E2F8E");
end
% cl
% linkaxes([ax1 ax2], 'x');
% xlim([1 800]);
% grid on;
% grid minor;
% 
% % (END) Symbol sync error function validation

