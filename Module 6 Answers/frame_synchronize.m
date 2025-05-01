function [message_bits, xc, peak_xcorr] = frame_synchronize(symbs, key, num_chars)
    % Compute cross-correlation with key
    xc = xcorr(symbs, key);
    xc = xc(length(symbs):end); % since the xcorr has length 2N - 1
    % xcorr pads the shorter vector with zeros at the front

    % Detect peak of xcorr
    n = length(xc);
    maxSum = -inf; % Initialize maxSum to negative infinity
    index = -1; % Initialize index to -1 in case no such index is found
    for i = 1:n

        if abs(xc(i)) > maxSum
            maxSum = abs(xc(i));
            index = i;
        end
    end
    peak_xcorr = index;
    message_bits = symbs(peak_xcorr + length(key):peak_xcorr + length(key) + num_chars * 8 - 1);

end