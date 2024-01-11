#include "nvwa/pctimer.h"
#include <iostream>

double pi(int n) {
    double sum = 0.0;
    int sign = 1;
    for (int i = 0; i < n; ++i) {           
        sum += sign/(2.0*i+1.0);
        sign *= -1;
    }
    return 4.0*sum;
}

int main(){
    auto t1 = nvwa::pctimer();
    double calculatedPi = pi(10000000);
    std::cout<<"Calculated PI = " << calculatedPi << std::endl;
    auto t2 = nvwa::pctimer();
    std::cout<<"Test took " << (t2 - t1) << " seconds"<<std::endl;
}