classdef CRC
    methods (Static)
        function remainder = remainder(input_bits)
            input_len = MAC.header_bits + MAC.msdu_bits;
            input = bit2int(input_bits', input_len);
            divisor_bits = [1 0 1 1];
            divisor_degree = MAC.crc_bits;
            divisor = bit2int(divisor_bits', divisor_degree + 1);
            divisor_shift = bitshift(divisor, (input_len - divisor_degree - 1) + divisor_degree);
           
            remainder = bitshift(input, divisor_degree);
            for k = 1:input_len
                if bitget(remainder, input_len + divisor_degree)
                    remainder = bitxor(remainder, divisor_shift);
                end
                remainder = bitshift(remainder, 1);
            end
            remainder = bitshift(remainder, -input_len);
        end

        function correct = check(input_bits, crc_remainder_bits)
            correct = false;
            if length(input_bits) == MAC.header_bits + MAC.msdu_bits
                remainder = CRC.remainder(input_bits);
                crc_remainder = bit2int(crc_remainder_bits', MAC.crc_bits);
                if remainder == crc_remainder
                    %disp('Message is error free.')
                    correct = true;
                end
            end
        end
    end
end
