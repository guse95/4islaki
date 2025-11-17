#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>


double wi(const int ind, const int n, const double* x, const double x0) {
    double res = 1;
    for (int i = 0; i < n; i++) {
        if (ind == i) continue;;
        res *= x0 - x[i];
    }
    return res;
}

double f(const double x) {
    // return std::log(x);
    // return std::sin(std::acos(-1) * x / 6);
    return std::asin(x);
}

double L(const double x, const int n, const double* xi) {
    double w_f[n] = {0}; // fi/wi
    double res = 0;
    for (int i = 0; i < n; i++) {
        w_f[i] = f(xi[i]) / wi(i, n, xi, xi[i]);

        res += w_f[i] * wi(i, n, xi, x);
    }
    return res;
}


std::vector<double> coef_Newton(const int n, const double* x) {
    std::vector table(n + 1, std::vector<double>(n + 1));

    for (int i = 0; i <= n; i++) {
        table[i][0] = f(x[i]);
    }

    for (int j = 1; j <= n; j++) {
        for (int i = 0; i <= n - j; i++) {
            table[i][j] = (table[i + 1][j - 1] - table[i][j - 1]) /
                         (x[i + j] - x[i]);
        }
    }

    std::vector<double> res(n + 1);
    for (int i = 0; i <= n; i++) {
        res[i] = table[0][i];
    }
    return res;
}


double P(double x, const int n, const double* xi) {
    auto coefs = coef_Newton(n - 1, xi);
    double res = coefs[0];
    double tmp = 1;
    for (int i = 1; i < n; i++) {
        tmp *= x - xi[i - 1];
        res += tmp * coefs[i];
    }
    return res;
}

int main() {

    int n = 4;
    // double xi[] = {0.1, 0.5, 0.9, 1.3};
    // double x_in = 0.8;

    // double xi[] = {0, 1, 2, 3};
    // double x_in = 1.5;

    // double xi[] = {-0.4, -0.1, 0.2, 0.5};
    double xi[] = {-0.4, 0, 0.2, 0.5};
    double x_in = 0.1;

    double val = f(x_in);
    double val_Lagrange = L(x_in, n, xi);
    std::cout << std::fixed << std::setprecision(5);
    std::cout << "Function: " << val << std::endl;
    std::cout << "Lagrange: " << val_Lagrange << std::endl;
    std::cout << "The absolute deviation of the interpolation is: " << fabs(val_Lagrange - val) << std::endl;

    double val_Newton = P(x_in, n, xi);
    std::cout << "Newton: " << val_Newton << std::endl;
    std::cout << "The absolute deviation of the interpolation is: " << fabs(val_Newton - val) << std::endl;
}