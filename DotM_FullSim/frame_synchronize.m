function [aligned_symbs, xc, peak_xcorr] = frame_synchronize(symbs, key, header)
    % Compute cross-correlation with key
    xc = xcorr(symbs, key);
    xc = xc(length(symbs):end); % since the xcorr has length 2N - 1
    % xcorr pads the shorter vector with zeros at the front

    % Detect peak of xcorr
    n = length(xc);
    maxSum = -inf; % Initialize maxSum to negative infinity
    index = -1; % Initialize index to -1 in case no such index is found
    for i = 1:(n-length(header))
        currentSum = abs(xc(i) + xc(i+length(key)) + xc(i+2*length(key)));
        currentSum = currentSum - sum(abs(xc(i+1:i+length(key)-1))) - sum(abs(xc(i+length(key)+1:i+2*length(key)-1)));
        if currentSum > maxSum
            maxSum = currentSum;
            index = i;
        end
    end
    peak_xcorr = index;
    aligned_symbs = symbs(peak_xcorr:end);

end