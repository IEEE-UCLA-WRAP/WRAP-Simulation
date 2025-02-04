function [I, Q, theta, err] = costas_loop(signal, fr, Fs, t, Kp, Ki, Kd)
    % Array Setup
    theta = zeros(1, length(signal));
    I = zeros(1, length(signal));
    Q = zeros(1, length(signal));
    I_ = zeros(1, length(signal));
    Q_ = zeros(1, length(signal));
    err = zeros(1, length(signal));
    
    % Design low-pass filter
    M = 5;
    lp = designfilt('lowpassfir', 'FilterOrder', M, 'CutoffFrequency', 1.5e6, 'SampleRate', Fs).Coefficients;
    
    % Iterate from index M + 1 to the end of the array
    for k = (M + 1):length(signal)

        % Mix the signal with the sines and cosines
        I_(k) = signal(k) * 2 * cos(2 * pi * fr * t(k) + theta(k));
        Q_(k) = signal(k) * -2 * sin(2 * pi * fr * t(k) + theta(k));
     
        % Lowpass with the previous M samples
        I_lowpassed = conv(lp, I_(k-M:k));
        Q_lowpassed = conv(lp, Q_(k-M:k));
        
        I(k) = I_lowpassed(M + 1);
        Q(k) = Q_lowpassed(M + 1);
        
        % Calculate the error at this sample point
        err(k) = I(k) * Q(k);
        
        % Update theta for next sample point
        theta(k+1) = theta(k) + Kp*err(k) + Ki*sum(err) + Kd*(err(k)-err(k-1));
    
    end
end