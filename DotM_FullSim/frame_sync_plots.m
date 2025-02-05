function frame_sync_plots(xc, aligned_symbols, start_symbol)
    %% Frame Synchronization Plotting
    plot(xc(start_symbol:end))
    title("Cross-Correlation with Header Key")
    figure
    
    stem(aligned_symbols(start_symbol:end))
    title("Frame Sync Symbols")
    figure
end