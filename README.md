# Wireless RF Analog Project

WRAP is an advanced project focused on wireless digital communications and RF circuits. WRAP teaches students a variety of circuits topics used in wireless transmitter and receiver  design, including amplifiers, mixers, and oscillators. On the signal processing side, WRAP covers the fundamentals of digital communication, digital filtering, and other techniques used in real-world communication systems. As the year progresses, students will use this knowledge to design, build, and test a physical wireless communication system. 

This repository contains all of the MATLAB simulation code used to model the communication system.

## Demo
https://github.com/IEEE-UCLA-WRAP/WRAP-Simulation/assets/30915871/6963a51b-d834-4fc1-8fde-54f530e71d5a

## High-Level System Design
A microcontroller on the transmitter side will encode a bit string onto a 1MHz carrier signal and pass it to the transmitter hardware, which upconverts it from 1 to 27MHz. This upconverted signal is then transmitted across the channel (air) and downconverted back to 1MHz by the receiver board. On the receiver side, another microcontroller will then interpret this signal and decode the transmitted bits. For now, our system emphasizes robustness over data throughput, so we implement BPSK as our digital modulations scheme, but we may change this in future iterations.

We model each of the communication steps one at a time in our MATLAB simulations. Listed below are the steps associated with the system transmitter, channel, and receiver.

### Transmitter Steps
- Map binary bits to symbols (e.g. BPSK +1 and -1), prepend packet header for frame synchronization, and upsample to get time-shifted delta functions.
- Filter delta functions with root-raised cosine filter to obtain baseband representation.
- Modulate baseband representation onto carrier wave (i.e. multiply baseband signal by complex exponential).
- Transmit passband signal (the real part).
  
### Channel
- Passband signal travels through wireless channel (acquires noise, experiences frequency and phase shifts, undergoes channel band-limit).
- 
### Receiver
- Receive passband signal.
- Demodulate received carrier signal using Costas loop carrier recovery algorithm, obtaining channel-modified received baseband signal.
- Perform matched filtering by applying another root-raised cosine filter.
- Sample received signal using zero-crossing timing recovery algorithm and perform symbol detection, obtaining received symbols.
- Cross-correlate symbols with key to detect start of data packet.
- Extract data symbols and map them back to binary bits.
