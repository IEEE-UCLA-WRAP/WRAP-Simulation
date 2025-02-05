function [symbs, tau, err] = timing_recovery_pll(x, sps, Kp, Ki, Kd, method)
%  TIMING_RECOVERY_PLL 
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
        
        err(i) = timing_error_detector(x, tau, i, sps, method);
        tau(i+1) = tau(i) + Kp*err(i) + Ki*sum(err) + Kd*(err(i)-err(i-1));
    
        i = i + 1;
    end
end