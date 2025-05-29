classdef PMD
    properties (Constant)
        sync = [1 1 1 -1 -1 -1 1 -1 -1 1 -1]; % Barker code
        sync_len = 11; % Barker code length
        ppdu_len = 48;  % # of bits, 6 chars * 8 bits/char

        Fc = 1e6; % Carrier frequency in Hz
        Fs = 4e6; % Sampling rate
        Rs = 50e3; % Symbol rate in symbols/sec (baud)
        B = 350e3; % Bandwidth
        num_symbols = 450; % Number of symbols/bits to send

        beta = 0.5; % RRC rolloff factor
        span = 10; % number of symbols for length of filter impulse response

        sps = PMD.Fs / PMD.Rs;
        num_samples = PMD.num_symbols * PMD.sps; % Number of discrete samples
        t = linspace(0, PMD.num_symbols / PMD.Rs, PMD.num_samples); % Time variable
        f = linspace(-PMD.Fs/2,PMD.Fs/2, PMD.num_samples); % Frequency variable        
        ps_filter = rcosdesign(PMD.beta, PMD.span, PMD.sps, 'sqrt'); % Use an SRRC as the pulse shape filter 
    end

    methods (Static)
        function transmitted_signal = transmitter(ppdu)
            % Add barker for synchronization
            bits = generate_rand_bits(325);
            barker = (PMD.sync + 1) * 0.5;
            frame = horzcat(barker, ppdu);
            for i = 1:10
                bits = horzcat(bits, frame);
            end
            bits = bits(1:PMD.num_symbols);

            % Convert bits to BPSK Symbols (1's and -1's)
            symbols = bits;
            symbols(symbols == 0) = -1;

            % Use upsample to form a delta train
            deltas = upsample(symbols, PMD.sps);

            % Convolve the deltas with the rectangular window
            transmitted_baseband = conv(deltas, PMD.ps_filter, 'same'); 
        
            % Modulate Baseband Signal to Passband
            phase_offset = pi/8;
            transmitted_signal = modulate_carrier(transmitted_baseband, PMD.Fc, PMD.t, phase_offset);
        end

        function ppdu = receiver(transmitted_signal)
            % Transmit Through Wireless Channel
            snr = 30;
            received_signal = add_channel_impairments(transmitted_signal, PMD.B, PMD.Fc, PMD.Fs, snr);
            
            %% Demodulation
            % Normalize the signal for stability of PID tuning constants
            normalized_signal = normalize(received_signal);

            % PID Tuning Constants
            Kp = 0.1;
            Ki = 0.002;
            Kd = 0;

            % Costas Loop to extract I and Q components
            [I, Q, theta, err] = costas_loop(normalized_signal, PMD.Fc, PMD.Fs, PMD.t, Kp, Ki, Kd);

            %% Apply Receiver Impulse Response

            % Combine I and Q components to construct the received baseband signal (I + jQ)
            received_baseband = I + 1j * Q;

            % Convolve the received baseband signal with the receiver filter
            received_samples = conv(received_baseband, PMD.ps_filter, 'same'); 

            %% Sample and Detect Symbols

            % TODO 1.1: Apply SPS frequency offset at the receiversps
            received_sps = PMD.sps + 5; % can adjust

            % TODO 1.2: Apply SPS phase offset at the receiver
            received_samples = [zeros(1, 25) received_samples]; % 23 to 52

            % Normalize for stability for PID stability
            received_samples = normalize(real(received_samples));

            % TODO Section 2: Implement one (or multiple) of the timing error detectors to properly sample
            % Recreate the downsample function iteratively, then implement the TED
            % PID Tuning Constants
            Kp = 6.8; % observe how Kp causes faster settling w/o oscillation but increases PMD-noise effect
            Ki = 1.1; % observe how Ki smooths the oscillations before settling to a constant offset (but eventually goes unstable)
            Kd = 0.1; % observe how Kd increases high frequency noise in the lock
            [symbs, sps_deltas] = downsample_new(received_samples, received_sps, Kp, Ki, Kd);

            sps_deltas_len = length(sps_deltas);
            sps_guess = zeros(1, sps_deltas_len) + received_sps;
            sps_guess_len = length(sps_guess);
            for i = 2:sps_guess_len
                sps_guess(i) = sps_guess(i) + (sps_deltas(i) - sps_deltas(i - 1));
            end

            % TODO: You'll probably want to change this once you have the timing error detector
            received_symbols = symbs;

            % Apply Thresholding to detect the symbols (turning them back to 1's and -1's)
            received_symbols(received_symbols > 0) = 1;
            received_symbols(received_symbols <= 0) = -1;

            %% Map Symbols back to Bits

            % Frame Sync Algorithm

            % TODO 3.1: Cross correlate and extract where the message starts
            [c, lags] = xcorr(received_symbols, PMD.sync);
            cross_corr = c(length(received_symbols):end);
            N = length(cross_corr);

            max_sum = -inf;
            ppdu_start = -1;
            for i = 1:N
                if cross_corr(i) > max_sum
                    max_sum = cross_corr(i);
                    ppdu_start = i;
                end
            end
            ppdu_start = ppdu_start + PMD.sync_len; % key length

            %% Extract bits
            rec_len = length(received_symbols);
            ppdu_bits = received_symbols;
            if ppdu_start < rec_len && ppdu_start + PMD.ppdu_len - 1 <= rec_len
                ppdu_bits = received_symbols(ppdu_start:ppdu_start + PMD.ppdu_len - 1);
            end
            for i = 1:length(ppdu_bits)
                if ppdu_bits(i) == -1
                    ppdu_bits(i) = 0;
                end
            end
            ppdu = char(strrep(strjoin(string(ppdu_bits)), " ", ""));
        end
    end
end