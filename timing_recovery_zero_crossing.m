function [recovered_signal, timing_error] = timing_recovery_zero_crossing(received_signal, sampling_rate)
    % received_signal: received signal samples
    % sampling_rate: sampling rate of the received signal

    % Compute the sampling period
    T = 1 / sampling_rate;

    % Compute the sign of the received signal samples
    sign_signal = sign(received_signal);

    % Initialize variables
    N = length(received_signal);
    recovered_signal = zeros(1, N);
    timing_error = zeros(1, N);

    % Loop through the samples
    for n = 2:N
        % Detect zero-crossings
        if sign_signal(n) ~= sign_signal(n - 1)
            % Estimate the timing error as the midpoint between the zero-crossings
            timing_error(n) = (n - 1) * T + T / 2;
        end
    end

    % Interpolate the timing error to estimate the timing offset for each sample
    interpolated_timing_error = interp1(find(timing_error), timing_error(timing_error ~= 0), 1:N, 'linear', 'extrap');

    % Adjust the sampling instants using the timing error
    for n = 1:N
        % Compute the index of the adjusted sampling instant
        adjusted_index = round(n - interpolated_timing_error(n) / T);

        % Ensure the adjusted index is within bounds
        adjusted_index = max(1, min(N, adjusted_index));

        % Assign the received signal sample to the adjusted sampling instant
        recovered_signal(n) = received_signal(adjusted_index);
    end
end
