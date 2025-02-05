function timing_offset_plots(transmit_bb, stored_transmit_bb, sps, m, timing_offset_plot_amount)
    %% Timing Offset Plots
    plot(transmit_bb)
    title("Transmit Baseband with Timing Offset")
    figure
    % The timing offset is minimally noticeable unless you zoom
    % in on a few symbols to compare to the previous Transmit
    % Baseband plot
    
    plot(stored_transmit_bb(length(m)*sps:length(m)*sps+timing_offset_plot_amount))
    title("Symbol-level View of Fractional Timing Offset")
    hold on
    plot(transmit_bb(length(m)*sps:length(m)*sps+1000))
    hold off
    legend("Without Offset", "With Offset")
    figure

    plot(stored_transmit_bb(length(m)*sps:length(m)*sps+timing_offset_plot_amount))
    title("Symbol-level View of Fractional Timing Offset")
    hold on
    plot(transmit_bb(length(m)*sps:length(m)*sps+timing_offset_plot_amount))
    
    stem(upsample(stored_transmit_bb(length(m)*sps:sps:length(m)*sps+timing_offset_plot_amount), sps))
    stem(upsample(transmit_bb(length(m)*sps:sps:length(m)*sps+timing_offset_plot_amount), sps))
    
    xline(0:sps:timing_offset_plot_amount,'--')
    
    hold off
    legend("Without Offset", "With Offset", "'Good' samples at high amplitudes on sps lines", "'Bad' samples close to zero on sps lines", "Lines Every 'sps'")
    figure

end