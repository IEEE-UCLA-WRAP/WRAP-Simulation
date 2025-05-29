%% Simulation Setup
alice = transmitter_node('00000001');
bob = receiver_node('00000101');
%alice.get_id()
chan = channel();
chan = chan.addTransmitter(alice);
chan = chan.addReceiver(bob);
%chan.getNodes();
alice.send_message(bob, chan);

% Test the following
% Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum.

