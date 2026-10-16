#include "fft.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <iostream>
using namespace hw1;
Vector dft(const Vector& x, bool inverse) {
    Vector y(x.size());
    for (std::size_t k=0; k<x.size(); ++k)
        for (std::size_t j=0; j<x.size(); ++j)
            y[k] += x[j] * std::polar(1.0, (inverse ? 2 : -2) * pi *
                double(k)*double(j)/double(x.size())) / (inverse ? double(x.size()) : 1.0);
    return y;
}
void close(const Vector& a, const Vector& b, double tol=1e-10) {
    for (std::size_t i=0; i<a.size(); ++i)
        if (std::abs(a[i]-b[i]) > tol) throw std::runtime_error("Numerical mismatch");
}
int main() {
    for (auto n : {1,2,3,7,16,31,32,64,65,127,128,256}) {
        Vector x(n);
        for (int i=0; i<n; ++i) x[i] = {std::sin(0.71*i)+0.2*i, std::cos(0.37*i)};
        close(fft(x), dft(x,false)); close(ifft(x), dft(x,true));
        close(fft(ifft(x)),x); close(ifft(fft(x)),x);
        double a=0,b=0;
        for (auto z:x) a+=std::norm(z);
        for (auto z:fft(x)) b+=std::norm(z)/n;
        if (std::abs(a-b)>1e-10*std::max(1.0,a)) throw std::runtime_error("Parseval failed");
        Vector impulse(n); impulse[0]=1; close(fft(impulse),Vector(n,1));
        close(fft(Vector(n)),Vector(n));
        std::cout << "PASS N=" << n << " DFT/IFFT/round-trip/Parseval/impulse/zero\n";
    }
    for (bool inverse : {false,true}) {
        bool rejected=false;
        try { auto unused = inverse ? ifft({}) : fft({}); (void)unused; }
        catch (const std::invalid_argument&) { rejected=true; }
        if (!rejected) throw std::runtime_error("Empty input accepted");
    }
    std::cout << "All tests passed.\n";
}
