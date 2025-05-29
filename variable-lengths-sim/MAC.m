classdef MAC
    properties (Constant)
        address_bits = 8;
        seq_control_bits = 6;
        header_bits = MAC.address_bits * 2 + MAC.seq_control_bits;
        msdu_bits = 16;
        crc_bits = 3;
    end
    methods (Static)
        function transmitted_signal = send_msdu(address1, address2, seq_control, msdu)
            mpdu = MAC.create_mpdu(address1, address2, seq_control, msdu);
            transmitted_signal = PLCP.send_psdu(mpdu); % MPDU = PSDU
        end

        function mpdu = create_mpdu(address1, address2, seq_control, msdu)
            addresses = horzcat(address1=='1', address2=='1');
            seq_control = int2bit(seq_control, MAC.seq_control_bits)';
            crc_remainder = CRC.remainder(horzcat(addresses, seq_control, msdu));
            mpdu = horzcat(addresses, seq_control, msdu, int2bit(crc_remainder, MAC.crc_bits)');
        end

        function [address1, address2, seq_control, msdu, err] = extract_msdu(transmitted_signal)
            % Check reults
            [mpdu, plcp_err] = PLCP.extract_psdu(transmitted_signal);
            address1 = '';
            address2 = '';
            seq_control = '';
            msdu = '';
            err = "ERR";
            if plcp_err == "OK"
                crc_correct = CRC.check(mpdu(1:end-3)=='1', mpdu(end-2:end)=='1');
                if crc_correct 
                    address1 = mpdu(1:MAC.address_bits);
                    address2 = mpdu(MAC.address_bits + 1:2 * MAC.address_bits);
                    seq_control = mpdu(2 * MAC.address_bits + 1:MAC.header_bits);
                    msdu = mpdu(MAC.header_bits + 1:end-3);
                    err = "OK";
                end
            end
        end
    end
end