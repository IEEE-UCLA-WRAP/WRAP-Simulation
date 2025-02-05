function symbs = timing_recovery(x, sps, Kp, Ki, Kd, method)
    % Array Setup
    i = 2;
    err = zeros(1, ceil(length(x) / sps));
    symbs = zeros(1, ceil(length(x) / sps));
    tau = zeros(1, 1 + ceil(length(x) / sps));
    % Draw first symbol -- no locking affects first symb, just draw naively
    symbs(1) = x(round(1*sps));
    
    % Begin locking
    while round(i*sps + tau(i)) < length(x)
        symbs(i) = x(round(i*sps + tau(i)));
        
        switch method
            case 'Early-Late Gate'
                err(i) = timing_error_detector_early_late_gate(x, tau, i, sps);
            case 'Gardner'
                err(i) = timing_error_detector_gardner(x, tau, i, sps);
            case 'Mueller-Muller'
                err(i) = timing_error_detector_mueller_muller(x, tau, i, sps);
            otherwise
                error('Invalid timing recovery method. Please provide either "Gardner", "Mueller-Muller", or "Early-Late Gate".');
        end
        tau(i+1) = tau(i) + Kp*err(i) + Ki*sum(err) + Kd*(err(i)-err(i-1));
    
        i = i + 1;
    end
end

% Mueller-Muller
function err = timing_error_detector_mueller_muller(x, tau, i, sps)
    curr_samp = x(round(i*sps + tau(i)));
    prev_samp = x(round((i - 1)*sps + tau(i-1)));
    err = sign(prev_samp)*curr_samp - sign(curr_samp)*prev_samp;
end

% Gardner
function err = timing_error_detector_gardner(x, tau, i, sps)
    zc_samp = x(round((i - 1/2)*sps + tau(i)));
    curr_samp = x(round(i*sps + tau(i)));
    prev_samp = x(round((i - 1)*sps + tau(i-1)));
    
    err = zc_samp * (prev_samp - curr_samp);
end

% Early-Late Gate 
function err = timing_error_detector_early_late_gate(x, tau, i, sps)
    symb_delta = 0.25; % Param to adjust!
    samp_delta = sps * symb_delta;

    early_samp = x(round(i*sps + tau(i) + samp_delta));
    prompt_samp = x(round(i*sps + tau(i)));
    late_samp = x(round(i*sps + tau(i) - samp_delta));
    err = prompt_samp * (early_samp - late_samp);
end