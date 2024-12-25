function [received_signal] = add_channel_impairments(signal, B, Fc, Fs, snr)
    % TODO: Apply Channel Band-Limit
    s = tf('s');
    Hc = B*s/(s^2 + s*B + Fc^2);
    Hd = c2d(Hc, 1/Fs, 'tustin');
    [a, b] = tfdata(Hd);
    %signal = real(filter(a{:}, b{:}, signal));
    
    % TODO: Add AWGN
    received_signal = awgn(signal, snr, 'measured');
end