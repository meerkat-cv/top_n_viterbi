# Top N Viterbi for OCR Decoding
A C++ code, without any dependencies and with a Python 3 wrapper, to do Top N Viterbi.

This code was transcripted from the Python implementation [here](https://github.com/carthach/kBestViterbi). Thanks for the greate work! It was very simple to follow and debug, however Python is not the fastest language around. From my use cases the runtime went from 850ms to ~5ms.

## Important Notes
This code was adapted for our needs, and the observation input from Viterbi was omitted. To be more precise the observation on time t is t, O(t)=t. It is quite trivial to back another input called for observation and do the required changes.

## TODO
- Change probabilities multiplication to log sums. We are working with fp32 so precision may become a problem.
- Add Viterbi function with observations input.
- Add tests.
