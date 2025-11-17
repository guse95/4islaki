#include <iostream>
#include <stdexcept>
#include <cmath>
#include <string>
#include "src/Matrix.h"


matrix<double> x_next_iter(matrix<double> &x){
    matrix<double> res = x;
    double x_1 = x(0, 0), x_2 = x(1, 0);
    res(0, 0) = sqrt(9 - x_2*x_2);
    res(1, 0) = log(x_1 + 3);
    return res;
}

double norm_phi_derivative(matrix<double> &x){
    int n = x.getRows();
    matrix<double> res = matrix<double>(n, n);
    double x_1 = x(0, 0), x_2 = x(1, 0);
    res(0, 0) = 0;
    res(0, 1) = -(2*x_2) / sqrt(9 - x_2*x_2);
    res(1, 0) = 1 / (x_1 + 3);
    res(1, 1) = 0;
    return std::max(res(1, 0), res(0, 1));
}

matrix<double> x_next_Newton(matrix<double> &x){
    matrix<double> res = x;
    matrix<double> J(2, 2), f(2, 1);
    double x_1 = x(0, 0), x_2 = x(1, 0);
    J(0, 0) = 2*x_1;
    J(0, 1) = 2*x_2;
    J(1, 0) = 1;
    J(1, 1) = -exp(x_2);
    matrix<double> J_inv = J.inverse();
    f(0, 0) = x_1*x_1 + x_2*x_2 - 9;
    f(1, 0) = x_1 - exp(x_2) + 3;
    matrix<double> J_inv_f = J_inv * f;
    res -= J_inv_f;
    return res;
}

void method_iter(matrix<double> (*x_next)(matrix<double>& ), double (*norm)(matrix<double>& ),double eps = 0.00001, int max_iter = 1000){
    matrix<double> x_prev(2, 1);
    x_prev(0,0) = 2.5;
    x_prev(1,0) = 1.5;

    matrix<double> x = x_prev;

    double q_l = fabs(norm(x));
    double q_r = fabs(norm(x_prev));
    double q;
    if (q_l > q_r) {
        q = q_l;
        x_prev = x;
    } else {
        q = q_r;
    }

    if (q >= 1.0){
        throw std::invalid_argument("The convergence condition is not satisfied norm = " + std::to_string(norm(x_prev)));
    }
    for (int i = 1; i <= max_iter; i++){
        x = x_next(x);
        matrix<double> diff = x - x_prev;
        if (std::max(fabs((diff(0, 0))), fabs(diff(1, 0))) * q / (1- q) < eps){
            std::cout << "Iteration method:" << "\n";
            std::cout << "X:\n" << x << "\n";
            std::cout << "Number of iterations: " << i << "\n";
            return;
        }
        x_prev = x;
    }
    throw std::runtime_error("Method iter did not converge in " + std::to_string(max_iter) + " iterations");
}

void method_Newton(matrix<double> (*x_next)(matrix<double>& ), double eps = 0.00001, int max_iter = 1000){

    matrix<double> x_prev(2, 1);
    x_prev(0,0) = 2.5;
    x_prev(1,0) = 1.5;

    matrix<double> x = x_prev;

    for (int i = 1; i <= max_iter; i++){
        x = x_next(x_prev);
        matrix<double> diff = x - x_prev;
        if (fabs(std::max(diff(0, 0), diff(1, 0))) < eps){
            std::cout << "Newton's method:" << std::endl;
            std::cout << "X:\n" << x << std::endl;
            std::cout << "Number of iterations: " << i << std::endl;
            return;
        }
        x_prev = x;
    }
    throw std::runtime_error("Method Newton did not converge in " + std::to_string(max_iter) + " iterations");
}

int main(){
    method_iter(&x_next_iter, &norm_phi_derivative);
    method_Newton(&x_next_Newton);
    return 0;
}