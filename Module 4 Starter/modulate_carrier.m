% TODO 1.1: Edit modulate_carrier to implement a constant phase_offset

function signal = modulate_carrier(transmited_baseband, Fc, t)

    carrier = exp(1j * (2 * pi * Fc * t));
    analytic = transmited_baseband .* carrier;
    signal = real(analytic);
end