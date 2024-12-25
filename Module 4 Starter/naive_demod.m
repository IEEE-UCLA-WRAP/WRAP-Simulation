function [I, Q] = naive_demod(signal, Fc, Fs, t)

    I_ = 2 * signal .* cos(2*pi*Fc*t);
    Q_ = 2 * signal .* -sin(2*pi*Fc*t);

    % Create a 1.5 MHz lowpass filter
    lp = designfilt('lowpassfir', 'FilterOrder', 6, 'CutoffFrequency', 1.5e6, 'SampleRate', Fs).Coefficients;

    % Analyze the generated filter
    filterAnalyzer(lp);

    % Lowpass filter the previous results to get I and Q
    I = conv(I_, lp, 'same');
    Q = conv(Q_, lp, 'same');

end