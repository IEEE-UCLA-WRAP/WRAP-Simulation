function costas_loop_plots(I, Q, theta, err, fc, frequency_variation, sps, m, phase_offset, start_index_plot, end_index_plot)
    %% Costas Loop Plotting
    plot(I(start_index_plot:end_index_plot))
    title("I")
    figure
    
    plot(Q(start_index_plot:end_index_plot))
    title("Q")
    figure
    
    plot(err(start_index_plot:end_index_plot))
    title("Err")
    figure
    
    plot(theta)
    title("Theta")
    figure
    
    plot(theta)
    hold on
    if frequency_variation ~= 0
        expected_theta = [zeros(1, length(m)*sps) 0:-(pi/2)*frequency_variation/(fc):(-(pi/2)*frequency_variation/(fc) * (length(theta)-length(m)*sps) - phase_offset)];
    else
        expected_theta = [zeros(1, length(m)*sps) phase_offset*ones(1, length(theta)-length(m)*sps)];
    end
    title("Theta vs Expected")
    plot(expected_theta)
    hold off
    figure
end