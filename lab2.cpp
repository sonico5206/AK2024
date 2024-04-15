#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <omp.h>
#include <string>

using namespace std;

//double f(double x) {
//    return x * x * x;
//}
//
//double calculateSecondDeriv(double x, double h) {
//    return (f(x - h) - 2 * f(x) + f(x + h)) / (h * h);
//}
//
//int main(int argc, char* argv[]) {
//    //int N = stoi(argv[1]);
//    int N = 5;
//    double a = -1.0;
//    double b = 3.0;
//    double h = (b - a) / N;
//    double sum = 0.0;
//
//    #pragma omp parallel for reduction(+: sum)
//    for (int i = 0; i <= N - 1; ++i) {
//        double x = a + i * h;
//        sum = max(sum, abs(calculateSecondDeriv(x, h)));
//    }
//
//    double accurateValue = 18.0;
//    double calculatedValue = sum;
//
//    cout << "Calculated value: " << calculatedValue << endl;
//    cout << "Accurate value: " << accurateValue << endl;
//    cout << "Difference: " << abs(calculatedValue - accurateValue) << endl;
//
//    return 0;
//}


double f(double x) {
    if (x >= 0 && x <= M_PI / 2) {
        return 1 - cos(x);
    }
    else {
        return x - M_PI / 2;
    }
}

int main(int argc, char* argv[]) {
    //int N = stoi(argv[1]);
    int N = 7;
    double a = 0.0;
    double b = 5.0;


    double accurateIntegral = 23/2 - M_PI*2 +  (M_PI * M_PI)/ 8; 
#pragma omp parallel for
        for (int i = 0; i < 4; i += 3) {
            int currentN = (i + 1) * N;

            double h = (b - a) / currentN;
            double sum = f(a) + f(b);
#pragma omp section
            {
#pragma omp parallel for reduction(+: sum)
                for (int j = 1; j < currentN; j += 2) {
                    double x = a + j * h;
                    sum += 4 * f(x);
                }
            }
#pragma omp section
            {
#pragma omp parallel for reduction(+: sum)
                for (int j = 2; j < currentN - 1; j += 2) {
                    double x = a + j * h;
                    sum += 2 * f(x);
                }
            }

            double calculatedIntegral = sum * h / 3;

            cout << "Calculated integral for N = " << currentN << ": " << calculatedIntegral << endl;
            cout << "Accurate integral: " << accurateIntegral << endl;
            cout << "Difference: " << abs(calculatedIntegral - accurateIntegral) << endl;
        }

    return 0;
}
