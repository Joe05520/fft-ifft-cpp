CXX = clang++
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra -Wpedantic
.PHONY: all test run clean
all: bin/ofdm
bin/ofdm: src/main.cpp src/fft.cpp src/fft.hpp
	mkdir -p bin
	$(CXX) $(CXXFLAGS) -Isrc src/main.cpp src/fft.cpp -o $@
bin/test_fft: tests/test.cpp src/fft.cpp src/fft.hpp
	mkdir -p bin
	$(CXX) $(CXXFLAGS) -Isrc tests/test.cpp src/fft.cpp -o $@
test: bin/test_fft
	./bin/test_fft
run: bin/ofdm
	./bin/ofdm all results
clean:
	rm -f bin/ofdm bin/test_fft
