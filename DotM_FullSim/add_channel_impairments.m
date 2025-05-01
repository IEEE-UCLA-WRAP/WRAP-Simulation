function [signal, phase_offset] = add_channel_impairments(signal, B, fc, Fs, noise_pwr)
    s = tf('s');
    Hc = B*s/(s^2 + s*B + fc^2);
    Hd = c2d(Hc, 1/Fs, 'tustin');
    [a, b] = tfdata(Hd);
    signal = real(filter(a{:}, b{:}, signal));
    
    % Non-ideal phase shift
    phase_offset = unifrnd(0, 2*pi);
    signal = real(signal*exp(1j*phase_offset));
    
    % Noise
    noise = noise_pwr * randn(1, length(signal));
    signal = signal + noise;

    plot(signal)
    title("Corrupted Signal")
    figure
end