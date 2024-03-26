function [recovered_signal, timing_error] = timing_recovery_gardner(received_signal, T, sampling_period)
    % received_signal: received signal samples
    % T: symbol period
    % sampling_period: sampling period of the received signal

    % Compute the delay for the early and late samples
    delay = round(T / (2 * sampling_period));

    % Initialize variables
    N = length(received_signal);
    yk = zeros(1, N);
    e = zeros(1, N);
    timing_error = zeros(1, N);
    recovered_signal = zeros(1, N);

    % Loop through the samples
    for n = delay+1:N-delay
        % Compute the early, late, and punctual samples
        x_early = received_signal(n - delay);
        x_late = received_signal(n + delay);
        x_punctual = received_signal(n);

        % Compute the error signal
        e(n) = (x_late - x_early) * x_punctual;

        % Accumulate the error signal
        yk(n) = yk(n - 1) + e(n);

        % Compute the timing error
        timing_error(n) = yk(n) * sampling_period^2;

        % Adjust the sampling instant using the timing error
        adjusted_index = round(n - timing_error(n) / sampling_period);
        adjusted_index = max(1, min(N, adjusted_index));
        recovered_signal(n) = received_signal(adjusted_index);
    end

end
