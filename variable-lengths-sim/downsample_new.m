
function [symbols, sps_deltas] = downsample_new(rec_samples, rec_sps, Kp, Ki, Kd)
    N = length(rec_samples);
    samp_width = floor(rec_sps * 0.1);

    symbols = zeros(1, ceil(N / rec_sps));
    sps_deltas = zeros(1, ceil(N / rec_sps));
    err = zeros(1, ceil(N / rec_sps));
    err_integral = zeros(1, ceil(N / rec_sps));
    
    i = 1;
    while floor(1 + (i - 1)*rec_sps + sps_deltas(i)) - samp_width < 1
        i = i + 1;
    end
    
    while floor(1 + (i - 1)*rec_sps + sps_deltas(i)) + samp_width <= N
        cur_t = floor(1 + (i - 1)*rec_sps + sps_deltas(i));
        symbols(i) = rec_samples(cur_t);

        prompt_val = rec_samples(cur_t);
        early_val = rec_samples(cur_t - samp_width);
        late_val = rec_samples(cur_t + samp_width);
        err(i) = sign(prompt_val) * (late_val - early_val);

        err_integral(i) = err_integral(i - 1) + err(i);
        err_deriv_i = err(i) - err(i - 1);
        sps_deltas(i + 1) = sps_deltas(i) + Kp*err(i) + Ki*err_integral(i) + Kd*err_deriv_i;

        i = i + 1;
    end

end