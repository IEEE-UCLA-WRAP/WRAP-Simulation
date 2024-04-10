# WRAP-Simulation

Matlab simulation of IEEE@UCLA WRAP DSP system, implementing BPSK modulation.

## High-level process
### Transmitter
- Map binary bits to symbols (e.g. BPSK +1 and -1), prepend packet header for frame synchronization, and upsample to get time-shifted delta functions.
- Filter delta functions with root-raised cosine filter to obtain baseband representation.
- Modulate baseband representation onto carrier wave (i.e. multiply baseband signal by complex exponential).
- Transmit passband signal (the real part).
### Channel
- Passband signal travels through wireless channel (acquires noise, experiences frequency and phase shifts, undergoes channel band-limit).
### Receiver
- Receive passband signal.
- Demodulate received carrier signal using Costas loop carrier recovery algorithm, obtaining channel-modified received baseband signal.
- Perform matched filtering by applying another root-raised cosine filter.
- Sample received signal using zero-crossing timing recovery algorithm and perform symbol detection, obtaining received symbols.
- Cross-correlate symbols with key to detect start of data packet.
- Extract data symbols and map them back to binary bits.

## Contributing

Create a new branch, titled with whatever new feature you're trying to implement.

## License

[MIT](https://choosealicense.com/licenses/mit/)
