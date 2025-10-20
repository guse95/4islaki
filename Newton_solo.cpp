#include <iostream>
#include <stdexcept>
#include <cmath>
#include <string>


double phi(double x){
    return sqrt((1.0 + std::log(x + 1)) / 2);
}

double phi_derivative(double x) {
    return 1 / (2 * sqrt(2) * (x + 1) * sqrt(1.0 + std::log(x + 1)));
}

double f(double x){
    return std::log(x + 1) - 2*x*x + 1;
}

double f_derivative(double x){
    return 1 / (x + 1) - 4 * x;
}

void method_iter(double (*phi)(double ), double (*phi_derivative)(double ),double eps = 0.00001, int max_iter = 1000){
    double x_prev = 1.5, x = 0;
    if (fabs(phi_derivative(x_prev)) >= 1.0){
        throw std::invalid_argument("The convergence condition is not satisfied");
    }
    for (int i = 1; i <= max_iter; i++){
        x = phi(x);
        if (fabs(x - x_prev) < eps){
            std::cout << "Iteration method:" << std::endl;
            std::cout << "X = " << x << std::endl;
            std::cout << "Last iter = " << i << std::endl << std::endl;
            return;
        }
        x_prev = x;
    }
    throw std::runtime_error("Method iter did not converge in " + std::to_string(max_iter) + " iterations");
}

void method_Newton(double (*f)(double ),double (*f_derivative)(double ) , double eps = 0.00001, int max_iter = 1000){
    double x_prev = 1.5, x = 0;

    for (int i = 1; i <= max_iter; i++){
        x = x_prev - f(x_prev) / f_derivative(x_prev);
        if (fabs(x - x_prev) < eps){
            std::cout << "Newton's method:" << std::endl;
            std::cout << "X = " << x << std::endl;
            std::cout << "Last iter = " << i << std::endl;
            return;
        }
        x_prev = x;
    }
    throw std::runtime_error("Method Newton did not converge in " + std::to_string(max_iter) + " iterations");
}

int main(){
    method_iter(&phi, &phi_derivative);
    method_Newton(&f, &f_derivative);
    return 0;
}