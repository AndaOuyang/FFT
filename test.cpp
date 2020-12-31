#include "fft.hpp"
#include <iostream>
int main(){
    std::vector<int> source_vec{1,5,-9,8,-7,5,6,-2,6,5,4,-2,-3,5,8,-9,-7,1,2,4};
    std::vector<std::complex<double>> vec(source_vec.begin(), source_vec.end());
    std::vector<std::complex<double>> freq = fft::fft(vec);
    std::vector<std::complex<double>> time = fft::ifft(freq);

    std::vector<int> round_to_int = fft::round(time);
    round_to_int.resize(source_vec.size());
    std::cout << ((round_to_int == source_vec) ? "fft then ifft recovers the time domain signal" : "fft failed") << '\n';
    return 0;
}