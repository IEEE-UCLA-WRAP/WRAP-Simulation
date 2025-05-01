function err = timing_error_detector_early_late_gate(x, tau, i, sps)
%  TIMING_ERROR_DETECTOR_EARLY_LATE_GATE 
        
    % Early-Late Gate Symbol-Fractional Delta
    symb_delta = 0.25; % Param to adjust!
    samp_delta = sps * symb_delta;
    % Early-Late Gate
    early_samp = x(round(i*sps + tau(i) + samp_delta));
    prompt_samp = x(round(i*sps + tau(i)));
    late_samp = x(round(i*sps + tau(i) - samp_delta));
    err = prompt_samp * (early_samp - late_samp);
end