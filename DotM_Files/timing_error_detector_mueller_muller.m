function err = timing_error_detector_mueller_muller(x, tau, i, sps)
%  TIMING_ERROR_DETECTOR_MUELLER_MULLER 
    % Mueller-Muller
    curr_samp = x(round(i*sps + tau(i)));
    prev_samp = x(round((i - 1)*sps + tau(i-1)));
    err = sign(prev_samp)*curr_samp - sign(curr_samp)*prev_samp;
end