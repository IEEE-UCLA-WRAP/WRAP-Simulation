% Clear all previous figures before starting
close all;

bits = [0 1 0 0 0 1 0 1 1 1 0 0 1 0 0 1 1 0 0 1 0 0 0 1 1 1 0 1 1 1 0 1 1 0,...
        0 0 0 1 0 1 1 1 0 0 1 0 0 1 1 0 0 1 0 1 0 0 1 0 0 0 0 0 0 1 0 1 0 0,...
        1 1 0 1 1 1 0 1 0 1 0 1 1 0 0 0 1 1 0 1 1 0 1 0 1 1 0 1 1 1 0 0 1 1,...
        0 0 1 0 0 0 0 0 0 1 0 1 0 0 1 1 0 1 1 1 0 1 1 1 0 1 1 0 1 0 0 1 0 1,...
        1 1 0 1 0 0 0 1 1 0 0 0 1 1 0 1 1 0 1 0 0 0 0 0 0 0 1 1 0 1 1 0 1 1,...
        1 1 0 1 1 0 1 1 0 1 0 1 1 1 0 0 0 0 0 1 1 1 0 1 0 1 0 1 1 1 0 1 0 0,...
        0 1 1 0 0 1 0 1 0 1 1 1 0 0 1 0 0 0 1 0 0 0 0 0 0 1 0 0 0 1 0 1 0 1,...
        1 0 1 1 1 0 0 1 1 0 0 1 1 1 0 1 1 0 1 0 0 1 0 1 1 0 1 1 1 0 0 1 1 0,...
        0 1 0 1 0 1 1 0 0 1 0 1 0 1 1 1 0 0 1 0 0 1 1 0 1 0 0 1 0 1 1 0 1 1,...
        1 0 1 0 1 0 0 1 0 1 1 0 1 0 1 1 1 1 1 1 1 1 0 0 0 1 0 0 1 0 0 1 0 1,...
        0 0 1 1 0 1 1 0 1 0 1 1 0 1 1 0 1 0 0 1 0 1 1 0 0 0 1 0 0 1 1 0 1 0,... 
        0 1 0 1 1 0 0 1 0 0 0 1 1 0 1 0 0 1 0 0 1 1 0 1 1 0 1 0 1 1 0 1 1 0,... 
        1 0 0 1 0 1 1 0 0 0 1 0 0 1 1 0 1 0 0 1 0 1 1 0 0 1 0 0 0 1 1 0 1 0,...
        1 0 0 1 1 0 0 1 1 1 0 0 1 0 0 0 0 1];

symbols = bits;
symbols(symbols == 0) = -1;

key = [+1 +1 +1 -1 -1 -1 +1 -1 -1 +1 -1];

% 01001000 01100001 

% Frame Sync Algorithm
num_chars = 7;
[message_bits, xc, peak_xcorr] = frame_synchronize(symbols, key, num_chars);
% Frame Sync Plots
figure;
plot(xc)
title("Cross-Correlation with Key")

%% Extract Binary message

% Map aligned symbols back to bits
message = message_bits;
message(message_bits == -1) = 0;

%% Decode and Display Message

s = string(message);
s = strjoin(s);
s = strrep(s, " ", "");
disp(s);
s_len = length(char(s));
inputString = char(s);
binaryString = inputString(1:end-mod(length(inputString),8));
binaryChunks = reshape(binaryString, 8, []).';
asciiChars = char(bin2dec(binaryChunks)).';
disp(asciiChars);
%% 
% Message should be "Congrats! Simulation is complete"
