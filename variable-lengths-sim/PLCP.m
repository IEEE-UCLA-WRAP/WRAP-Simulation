classdef PLCP
    properties (Constant)
        barker = [1 1 1 -1 -1 -1 1 -1 -1 1 -1];
        len_bits = 6;
    end
    methods (Static)
        function ppdu = create_ppdu(psdu)
            psdu_len = dec2bin(length(psdu), PLCP.len_bits)=='1';
            parity_bit = dec2bin(mod(nnz(psdu_len), 2))=='1';
            ppdu = horzcat(psdu_len, parity_bit, psdu);
        end

        function transmitted_signal = send_psdu(psdu)
            ppdu = PLCP.create_ppdu(psdu);
            transmitted_signal = PMD.transmitter(ppdu);
        end

        function [mpdu, err] = extract_psdu(transmitted_signal)
            ppdu = PMD.receiver(transmitted_signal);
            psdu_len = ppdu(1:PLCP.len_bits);
            parity_bit = ppdu(7);
            mpdu = '';
            err = "ERR";
            if dec2bin(mod(nnz(psdu_len=='1'), 2)) == parity_bit
                mpdu = ppdu(8:end);
                err = "OK";
            end 
        end        
    end
end

% frame_payloads = [[0 1 0 0 0 1 0 1 1 1 0 0 1 0 0 1 1 0 0 1 0 0 0 1 1 1 0 1 1 1 0 1 1 0,...
%                     0 0 0 1 0 1 1 1 0 0 1 0 0 1 1 0 0 1 0 1 0 0 1 0 0 0 0 0 0 1 0 1 0,... 
%                     0 1 1 0 1 1 1 0,... % n
%                     1 0 1 0 1 1 0 0,... %
%                     0 1 1 0 1 1 0 1,... % m
%                     0 1 1 0 1 1 1 0,... % n
%                     0 1 1 0 0 1 0 0 0 0 0,... 
%                     0 1 0 1 0 0 1 1,... % S
%                     0 1 1 1 0 1 1 1,... % w
%                     0 1 1 0 1 0 0 1,... % i
%                     0 1 1 1 0 1 0 0,... % t
%                     0 1 1 0 0 0 1 1,... % c
%                     0 1 1 0 1 0 0 0,... % h
%                     0 0 0 0 1 1 0 1 1 0 1 1,...
%                     1 1 0 1 1 0 1 1 0 1 0 1 1 1 0 0 0 0 0 1 1 1 0 1 0 1 0 1 1 1 0 1 0 0 0,...
%                     1 1 0 0 1 0 1 0,...
%                     1 1 1 0 0 1 0,... % 7-bit Barker code
%                     0 0 1 0 0 0 0 0,... %
%                     0 1 0 0 0 1 0 1,... % E
%                     0 1 1 0 1 1 1 0,... % n
%                     0 1 1 0 0 1 1 1,... % g
%                     0 1 1 0 1 0 0 1,... % i
%                     0 1 1 0 1 1 1 0,... % n
%                     0 1 1 0 0 1 0 1,... % e
%                     0 1 1 0 0 1 0 1,... % e
%                     0 1 1 1 0 0 1 0,... % r
%                     0 1 1 0 1 0 0 1,... % i
%                     0 1 1 0 1 1 1 0,... % n
%                     1 0 1 0 0 1 0 1 1 0 1 0 1 1 1 1 1,...
%                     1 1 1 0 0 0 1 0 0 1 0,... % 11-bit barker code
%                     0 1 0 1 0 0 1 1,... % S
%                     0 1 1 0 1 0 1 1,... % k
%                     0 1 1 0 1 0 0 1,... % i
%                     0 1 1 0 0 0 1 0,... % b
%                     0 1 1 0 1 0 0 1,... % i
%                     0 1 1 0 0 1 0 0,... % d
%                     0 1 1 0 1 0 0 1,... % i
%                     0 0 1 1 0 1 1 0 1 0 1 1 0 1 1 0 1 0 0 1 0 1 1 0 0 0 1 0 0 1 1 0 1 0,...
%                     0 1 0 1 1 0 0 1 0 0 0 1 1 0 1 0 1 0 0 1 1 0 0 1 1 1 0 0 1 0 0 0 0 1,...
%                 ]];
