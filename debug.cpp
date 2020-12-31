#include <complex>
#include <iostream>

int main(){
    std::complex<double> c1 = {1, 5.3};
    c1 *= 5;
    std::cout << c1 << '\n';
    return 0;
}