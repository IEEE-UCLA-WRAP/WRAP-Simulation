function [I, Q, theta, err] = costas_loop(signal, fr, Fs, t, Kp, Ki, Kd)
    % Array Setup
    theta = zeros(1, length(signal));
    I = zeros(1, length(signal));
    Q = zeros(1, length(signal));
    I_ = zeros(1, length(signal));
    Q_ = zeros(1, length(signal));
    err = zeros(1, length(signal));
    
    % TODO 3.1.0: Design low-pass filter
    M = 5;
    lp = designfilt('lowpassfir', 'FilterOrder', M, 'CutoffFrequency', 1.5e6, 'SampleRate', Fs).Coefficients;
    
    % Iterate from index M + 1 to the end of the array
    for k = (M + 1):length(signal)

        % TODO 3.1.1: Mix the signal with the sines and cosines
        I_(k) = signal(k) * 2 * cos(2 * pi * fr * t(k) + theta(k));
        Q_(k) = signal(k) * -2 * sin(2 * pi * fr * t(k) + theta(k));
     
        % TODO 3.1.2: Lowpass the past M samples
        % Hint: use conv with the "full" parameter, then extract the M+1th
        % element from the resulting array
        I_lowpassed = conv(lp, I_(k-M:k));
        Q_lowpassed = conv(lp, Q_(k-M:k));
        
        I(k) = I_lowpassed(M + 1);
        Q(k) = Q_lowpassed(M + 1);
        
        % TODO 3.2.1: Calculate the error at this sample point
        err(k) = I(k) * Q(k);
        
        % TODO 3.2.2: Update theta at this sample point
        theta(k+1) = theta(k) + Kp*err(k) + Ki*sum(err) + Kd*(err(k)-err(k-1));
    
    end
end