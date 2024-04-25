function [I, Q, theta, err] = carrier_recovery_pll(signal, fr, Fs, t, Kp, Ki, Kd)
%  CARRIER_RECOVERY_PLL 
    % Array Setup
    theta = zeros(1, length(signal));
    I = zeros(1, length(signal));
    Q = zeros(1, length(signal));
    I_mod = zeros(1, length(signal));
    Q_mod = zeros(1, length(signal));
    err = zeros(1, length(signal));
    
    % Design low-pass filter
    order = 5;
    lp = designfilt('lowpassfir', 'FilterOrder', order, 'CutoffFrequency', 0.2 * fr / Fs);
    
    for i = (1 + order):length(signal)
        % Costas Loop Joint Demodulation & Carrier Phase Recovery
        %1) Demodulate
        I_mod(i) = signal(i) * 2 * cos(2 * pi * fr * t(i) + theta(i));
        Q_mod(i) = signal(i) * -2 * sin(2 * pi * fr * t(i) + theta(i));
     
        %2) Lowpass
        Ismoothed = filter(lp, I_mod(i-order:i));
        Qsmoothed = filter(lp, Q_mod(i-order:i));
        
        %3) Extract I and Q
        I(i) = Ismoothed(end);
        Q(i) = Qsmoothed(end);
        
        %4) Calculate Error (Costas Phase Error Detector)
        err(i) = I(i) * Q(i);
        
        %5) Update Theta
        theta(i+1) = theta(i) + Kp*err(i) + Ki*sum(err) + Kd*(err(i)-err(i-1));
    
    end
end