function signal = modulate_carrier(transmitted_baseband, Fc, t, phase_offset)
    carrier = exp(1j * (2 * pi * Fc * t + phase_offset));
    analytic = transmitted_baseband .* carrier;
    signal = real(analytic);
end