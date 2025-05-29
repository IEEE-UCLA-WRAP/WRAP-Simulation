classdef transmitter_node
    properties
        my_id
        current_seq;
        channel
    end
    methods
        function node = transmitter_node(id)
            node.my_id = id;
            node.current_seq = 0;
        end

        function id = get_id(node)
            id = node.my_id;
        end

        function send_message(node, dest_node, channel)
            % Prompt for the message send
            prompt = "Input your message:";
            full_msg = inputdlg(prompt, "s");
            fprintf("Sending: %s\n", full_msg{1,1})
            full_msg_len = length(full_msg{1,1});

            % Split message into appropriate number of payloads
            payload_len = MAC.msdu_bits / 8;
            payload = zeros(1, MAC.msdu_bits);

            ch_count = payload_len;
            for i=1:full_msg_len
                ch = full_msg{1,1}(i);
                ch_arr = dec2bin(ch, 8)=='1'; % Binary of the ascii value.
                for j = 1:8
                    payload(8*(payload_len - ch_count)+j) = ch_arr(j);
                end
                ch_count = ch_count - 1;

                % Transmit payload 
                if ch_count == 0 || (i == full_msg_len && ch_count > 0)
                    address1 = node.my_id;
                    address2 = dest_node.get_id;
                    seq_control = node.current_seq;
                    for times = 1:4 % Send multiple times due to rather high possibility of error
                        transmitted_signal = MAC.send_msdu(address1, address2, seq_control, payload);
                        channel.broadcast(transmitted_signal);
                    end
                    if node.current_seq < 63
                        node.current_seq = node.current_seq + 1;
                    else
                        node.current_seq = 0;
                    end
                    
                    % reset
                    ch_count = payload_len;
                    payload = zeros(1,ch_count*8);
                end
            end
        end
    end
end