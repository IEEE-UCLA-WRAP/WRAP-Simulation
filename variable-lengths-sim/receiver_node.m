classdef receiver_node
    properties
        my_id
        current_seqs
        channel
    end
    methods
        function node = receiver_node(id)
            node.my_id = id;
            node.current_seqs = containers.Map('KeyType','char','ValueType','double');
        end

        function id = get_id(node)
            id = node.my_id;
        end

        function receive_message(node, transmitted_signal, channel)
            allowed_skips = 3;
            [address1, address2, seq_control, payload, mac_error] = MAC.extract_msdu(transmitted_signal);

            if mac_error == "OK"
                if address2 ~= node.my_id
                    return
                end

                seq = bin2dec(seq_control);
                if isKey(node.current_seqs, address1)
                    prev_seq = node.current_seqs(address1);
                    if (seq == 0) && (prev_seq > 0)
                        node.current_seqs(address1) = 0;
                    elseif (seq - prev_seq > 0) && (seq - prev_seq < allowed_skips)
                        node.current_seqs(address1) = seq;
                    elseif (seq + (63 - prev_seq)) < allowed_skips
                        node.current_seqs(address1) = seq;
                    else
                        return
                    end
                else
                    if channel.isTransmitter(address1) == true
                        node.current_seqs(address1) = 0;
                    else
                        return
                    end
                end
                
                s_len = length(char(payload));
                inputString = char(payload);
                binaryString = inputString(1:end-mod(s_len,8));
                binaryChunks = reshape(binaryString, 8, []).';
                asciichars = char(bin2dec(binaryChunks)).';
                %fprintf("%s : %s\n", payload, asciichars);
                fprintf("%s", asciichars);
            end
        end
    end
end