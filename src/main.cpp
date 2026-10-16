#include "fft.hpp"
#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace hw1;
Vector pattern(std::size_t n, int p) {
    Vector a(n);
    const double levels[] = {-3, -1, 1, 3}; // Unnormalized square 16-QAM.
    std::uint32_t state = 313512050u;
    auto next = [&]() { state = 1664525u * state + 1013904223u; return state; };
    for (std::size_t k = 0; k < n; ++k) {
        if (p == 0) a[k] = {1, 1};
        if (p == 1) a[k] = k % 2 ? Complex(-3, -1) : Complex(3, 1);
        if (p == 2) a[k] = {levels[k % 4], levels[(k / 4) % 4]};
        if (p == 3) {
            const auto re = (next() >> 30) & 3u;
            const auto im = (next() >> 30) & 3u;
            a[k] = {levels[re], levels[im]};
        }
    }
    return a;
}
double nearest(double x) {
    const double levels[] = {-3, -1, 1, 3};
    return *std::min_element(std::begin(levels), std::end(levels),
        [x](double a, double b) { return std::abs(x-a) < std::abs(x-b); });
}
void run(std::size_t n, const std::filesystem::path& out, std::ostream& summary) {
    for (int p = 0; p < 4; ++p) {
        const auto x = pattern(n, p);
        const auto spectrum = fft(x); // Standalone FFT benchmark.
        const auto tx = ifft(x);     // OFDM transmitter.
        const auto channel = tx;    // Ideal unit-gain, noiseless channel.
        const auto rx = fft(channel);
        double maxerr = 0, energy = 0, errenergy = 0;
        std::size_t errors = 0;
        std::ofstream f(out / ("N" + std::to_string(n) + "_pattern" + std::to_string(p+1) + ".csv"));
        if (!f) throw std::runtime_error("Cannot create result CSV");
        f << "index,input_re,input_im,fft_re,fft_im,ifft_re,ifft_im,rx_re,rx_im\n" << std::setprecision(17);
        for (std::size_t k = 0; k < n; ++k) {
            maxerr = std::max(maxerr, std::abs(rx[k]-x[k]));
            energy += std::norm(x[k]); errenergy += std::norm(rx[k]-x[k]);
            errors += Complex(nearest(rx[k].real()), nearest(rx[k].imag())) != x[k];
            f << k << ',' << x[k].real() << ',' << x[k].imag() << ','
              << spectrum[k].real() << ',' << spectrum[k].imag() << ','
              << tx[k].real() << ',' << tx[k].imag() << ','
              << rx[k].real() << ',' << rx[k].imag() << '\n';
        }
        summary << n << ',' << p+1 << ',' << maxerr << ',' << std::sqrt(errenergy/energy)
                << ',' << errors << ',' << n << '\n';
        if (maxerr > 1e-10 || errors) throw std::runtime_error("OFDM verification failed");
        std::cout << "N=" << n << " pattern=" << p+1 << " max_error=" << maxerr
                  << " symbol_errors=" << errors << '/' << n << '\n';
    }
}
int main(int argc, char** argv) {
    try {
        if (argc > 3) throw std::invalid_argument("Usage: ofdm [N|all] [output_directory]");
        const std::string arg = argc > 1 ? argv[1] : "all";
        std::vector<std::size_t> sizes = {3, 32, 64, 127, 128};
        if (arg != "all") {
            if (arg.empty() || arg.find_first_not_of("0123456789") != std::string::npos)
                throw std::invalid_argument("N must be a positive decimal integer");
            const auto n = std::stoull(arg);
            if (!n || n > 1048576) throw std::invalid_argument("Demo supports 1 <= N <= 1048576");
            sizes = {static_cast<std::size_t>(n)};
        }
        const std::filesystem::path out = argc > 2 ? argv[2] : "results";
        std::filesystem::create_directories(out);
        std::ofstream summary(out / "cpp_summary.csv");
        if (!summary) throw std::runtime_error("Cannot create summary");
        summary << "N,pattern,ofdm_max_error,evm_rms,symbol_errors,symbols\n" << std::setprecision(17);
        for (auto n : sizes) run(n, out, summary);
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
