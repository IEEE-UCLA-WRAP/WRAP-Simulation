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
    lp = ?
    
    % Iterate from index M + 1 to the end of the array
    for k = (M + 1):length(signal)

        % TODO 3.1.1: Mix the signal with the sines and cosines
        I_(k) = ?
        Q_(k) = ?
     
        % TODO 3.1.2: Lowpass the past M samples
        % Hint: use conv with the "full" parameter, then extract the M+1th
        % element from the resulting array
        
        I(k) = ? % (Does not have to be a one line solution)
        Q(k) = ?
        
        % TODO 3.2.1: Calculate the error at this sample point
        % err(k) = ?;
        
        % TODO 3.2.2: Update theta at this sample point
        % theta( ? ) = ?;
    
    end
end