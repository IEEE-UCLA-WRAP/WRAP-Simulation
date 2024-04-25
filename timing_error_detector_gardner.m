function err = timing_error_detector_gardner(x, tau, i, sps)
%  TIMING_ERROR_DETECTOR_GARDNER 
    % Gardner
    zc_samp = x(round((i - 1/2)*sps + tau(i)));
    curr_samp = x(round(i*sps + tau(i)));
    prev_samp = x(round((i - 1)*sps + tau(i-1)));
    
    err = zc_samp * (prev_samp - curr_samp);
end