# Validation record

Validated locally on 2026-10-16 using C++17 / Apple Clang and MATLAB R2024a.

- `make all test`: passed independent C++ direct DFT/IFFT, round-trip, Parseval, impulse, and zero-input checks at N = 1, 2, 3, 7, 16, 31, 32, 64, 65, 127, 128, 256; empty input rejected.
- Fresh C++ ideal-channel OFDM runs: four patterns at N = 3, 32, 64, 127, 128; zero symbol errors in every case.
- `verify_fft`: 20 retained cases passed independent MATLAB FFT/IFFT/reconstruction comparisons with a maximum absolute error threshold of 1e-10.
- `plot_fft_figures`: regenerated all 10 figures from retained measurements with MATLAB; exported PNGs were visually reviewed.
- Core C++ implementation and existing C++ tests are byte-identical to the source implementation used to prepare this public snapshot.

GitHub Actions independently builds and tests the C++ implementation on Linux and macOS. MATLAB verification is a separate local check; it is not part of CI.
