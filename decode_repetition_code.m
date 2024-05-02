function message_bits = decode_repetition_code(coded_bits, n)
    num_bits = length(coded_bits)/n;
    message_bits = mode(reshape(coded_bits, n, num_bits),1);
end