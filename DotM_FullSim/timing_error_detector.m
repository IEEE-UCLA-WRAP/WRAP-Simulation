function err = timing_error_detector(x, tau, i, sps, method)
    switch method
        case 'Early-Late Gate'
            err = timing_error_detector_early_late_gate(x, tau, i, sps);
        case 'Gardner'
            err = timing_error_detector_gardner(x, tau, i, sps);
        case 'Mueller-Muller'
            err = timing_error_detector_mueller_muller(x, tau, i, sps);
        otherwise
            error('Invalid timing recovery method. Please provide either "Gardner", "Mueller-Muller", or "Early-Late Gate".');
    end
end