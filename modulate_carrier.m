function signal = modulate_carrier(transmit_bb, fc, t)
    carrier = exp(1j*2*pi*fc * t);
    
    analytic = transmit_bb .* carrier(1:length(transmit_bb));
    
    signal = real(analytic);

    plot(signal)
    title("Signal")
    figure
end