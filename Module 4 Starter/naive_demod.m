function [I, Q] = naive_demod(signal, Fc, Fs, t)

    I_ = 2 * signal .* cos(2*pi*Fc*t);
    Q_ = 2 * signal .* -sin(2*pi*Fc*t);

    % TODO 2.1.0: Remove the below 2 lines so we can replace with our own lowpass filter
    I = lowpass(I_, 1.5e6, Fs);
    Q = lowpass(Q_, 1.5e6, Fs);

    % % TODO 2.1.1: Generate a 1.5 MHz lowpass filter
    % lp = ?

    % % TODO 2.1.2: Analyze the generated filter (just check if bode plot looks correct)
    % filterAnalyzer(lp);

    % % TODO 2.1.3: Convolve filter with previous results to get I and Q
    % I = ?
    % Q = ?

end