function timing_recovery_plots(timing_offset_samps, received_samples, x, sps, sps_offset, m, symbs, tau, err, start_symbol_plot, end_symbol_plot, start_symbol_constellation, end_symbol_constellation)
    %% Plotting for Timing Recovery
    sps = round(sps - sps_offset);
    
    plot(x(start_symbol_plot*sps:end_symbol_plot*sps))
    title("Input Signal")
    figure
    
    plot(err(start_symbol_plot:end_symbol_plot))
    title("Err")
    figure
    
    plot(symbs(start_symbol_plot:end_symbol_plot))
    title("Timing Recovered Symbols")
    figure
    
    stem(symbs(start_symbol_plot:end_symbol_plot))
    title("Timing Recovered Symbols")
    figure
    
    scatterplot(symbs(start_symbol_constellation:end_symbol_constellation))
    title("Timing Recovered Symbol Constellation")
    figure
    
    scatterplot(received_samples(start_symbol_constellation*sps:sps:end_symbol_constellation*sps))
    title("Non-Timed Symbol Constellation")
    figure
    
    plot(x(start_symbol_plot*sps:end_symbol_plot*sps))
    hold on
    interp_err = interp1(sps:sps:length(err)*sps, err, 0:length(err)*sps, 'linear');
    plot(interp_err(start_symbol_plot*sps:end_symbol_plot*sps))
    title("Overlay of Error")
    hold off
    figure
    
    plot(x(start_symbol_plot*sps+1:end_symbol_plot*sps))
    hold on
    interp_symbs = interp1(sps:sps:length(symbs)*sps, symbs, 0:length(symbs)*sps, 'linear');
    plot(interp_symbs((start_symbol_plot-1)*sps+1-timing_offset_samps:(end_symbol_plot-1)*sps-timing_offset_samps))
    title("Overlay of Symbols")
    legend("Samples", "Interpolated Symbols")
    hold off
    figure

    plot(tau)
    hold on
    if sps_offset ~= 0
        expected_tau = [zeros(1, length(m)) 0:-sps_offset:(-sps_offset * (length(tau)-length(m)) - timing_offset_samps)];
    else
        expected_tau = [zeros(1, length(m)) timing_offset_samps * ones(1, length(tau)-length(m))];
    end
    title("Tau vs Expected")
    plot(expected_tau)
    hold off
    figure
    
    plot(x(start_symbol_plot*sps:end_symbol_plot*sps))
    hold on
    x_downsample = [downsample(x,sps) zeros(1,30)];
    interp_raw = interp1(0:sps:(length(symbs)-1)*sps, x_downsample(1:length(symbs)), 0:(length(symbs))*sps, 'linear');
    plot(interp_raw(start_symbol_plot*sps-timing_offset_samps:end_symbol_plot*sps-timing_offset_samps))
    %x_up = upsample(x_downsample,sps);
    %stem(x_up(start_symbol_plot*sps:end_symbol_plot*sps))
    title("Overlay of Raw Downsamples")
    legend("Samples", "Interpolated Symbols")
    hold off
end