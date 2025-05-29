classdef channel
    properties
        transmitter_nodes = []
        receiver_nodes = []
    end

    methods 
        function channel = channel()
            channel.transmitter_nodes = [];
            channel.receiver_nodes = [];
        end

        function channel = addTransmitter(channel, node) 
            channel.transmitter_nodes = [channel.transmitter_nodes; node];
        end

        function channel = addReceiver(channel, node) 
            channel.receiver_nodes = [channel.receiver_nodes; node];
        end

        function [transmitter_nodes, receiver_nodes] = getNodes(channel) 
            transmitter_nodes = channel.transmitter_nodes;
            receiver_nodes = channel.receiver_nodes;
        end

        function status = isReceiver(channel, id) 
            status = false;
            n_count = length(channel.receiver_nodes);
            for i = 1:n_count
                if id == channel.receiver_nodes(i).get_id()
                    status = true;
                    return
                end
            end
        end

        function status = isTransmitter(channel, id) 
            status = false;
            n_count = length(channel.transmitter_nodes);
            for i = 1:n_count
                if id == channel.transmitter_nodes(i).get_id()
                    status = true;
                    return
                end
            end
        end

        function broadcast(channel, transmitted_signal)
            n_count = length(channel.receiver_nodes);
            for i = 1:n_count
                channel.receiver_nodes(i).receive_message(transmitted_signal, channel);
            end
        end
    end
end