#include <iostream>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include "../src/Matrix.h"

using vd = std::vector<double>;

vd run_throw(const int n, matrix<double>& matrixA, vd& d) {
    vd P(n);
    vd Q(n);
    vd X(n);
    P[0] = (-1) * matrixA[0][1]/matrixA[0][0];
    Q[0] = d[0] / matrixA[0][0];

    for (int i = 1; i < n - 1; ++i) {
        double denom = matrixA[i][i] + matrixA[i][i - 1] * P[i - 1];
        P[i] = (-1) * matrixA[i][i + 1] / denom;
        Q[i] = (d[i] - matrixA[i][i - 1] * Q[i - 1]) / denom;
    }

    double denom = matrixA[n - 1][n - 1] + matrixA[n - 1][n - 2] * P[n - 2];
    Q[n - 1] = (d[n - 1] - matrixA[n - 1][n - 2] * Q[n - 2]) / denom;

    X[n - 1] = Q[n - 1];
    for (int i = n - 2; i >= 0; i--) {
        X[i] = P[i] * X[i + 1] + Q[i];
    }
    return X;
}

void initializeCoefficients(const int n, const vd& y, const vd& h, const std::vector<double>& c_solution,
     vd& a, vd& b, vd& c, vd& d) {

    c[0] = 0.0;
    for (int i = 1; i <= n; i++) {
        if (i - 1 < c_solution.size()) {
            c[i] = c_solution[i - 1];
        } else {
            c[i] = 0.0;
        }
    }

    for (int i = 0; i <= n; i++) {
        a[i] = y[i];

        if (i < n) {
            b[i] = ((y[i + 1] - y[i]) / h[i]) -
                   (h[i] / 3.0) * (c[i + 1] + 2.0 * c[i]);
            d[i] = (c[i + 1] - c[i]) / (3.0 * h[i]);
        } else {
            b[i] = (y[i + 1] - y[i]) / h[i] -
                   (2.0 * h[i] / 3.0) * c[i];
            d[i] = -c[i] / (3.0 * h[i]);
        }
    }
}

void calculateSplineCoefficients(const int n, const vd& x, const vd& y, vd& a, vd& b, vd& c, vd& d) {
    //h_i = x_i - x_{i-1}
    vd h(n);

    for (int i = 1; i <= n; i++) {
        h[i - 1] = x[i] - x[i - 1];
    }

    matrix<double> A(n - 1, n - 1);
    vd rght(n - 1);

    if (n >= 2) {
        //2(h1 + h2)c2 + h2c3 = 3[(f2 - f1)/h2 - (f1 - f0)/h1]
        A[0][0] = 2.0 * (h[0] + h[1]);
        if (n > 2) {
            A[0][1] = h[1];
        }
        rght[0] = 3.0 * ((y[2] - y[1]) / h[1] -
                     (y[1] - y[0]) / h[0]);

        for (int i = 2; i <= n - 1; i++) {
            A[i - 1][i - 2] = h[i - 1]; // h_{i-1} для c_{i-1}

            A[i - 1][i - 1] = 2.0 * (h[i - 1] + h[i]); // 2(h_{i-1} + h_i) для c_i

            if (i < n - 1) {
                A[i - 1][i] = h[i]; // h_i для c_{i+1}
            }

            rght[i - 1] = 3.0 * ((y[i + 1] - y[i]) / h[i] - (y[i] - y[i - 1]) / h[i - 1]);
        }

        std::vector<double> c_solution = run_throw(n - 1, A, rght);

        //a, b, d
        initializeCoefficients(n, y, h, c_solution, a, b, c, d);
    } else {
        //n = 1
        initializeCoefficients(n, y, h, vd(n, 0.0), a, b, c, d);
    }
}

double evaluate(const int n, const double x, vd& xi, vd& a, vd& c, vd& b, vd& d) {
    // Находим интервал, которому принадлежит x
    int segment = 0;
    for (int i = 0; i <= n; i++) {
        if (x <= xi[i]) {
            segment = i - 1;
            break;
        }
    }
    if (segment > n) segment = n;

    double dx = x - xi[segment];
    return a[segment] + b[segment] * dx + c[segment] * dx * dx + d[segment] * dx * dx * dx;
}

int main() {
    const int n = 4;
    vd a(n);
    vd b(n);
    vd c(n);
    vd d(n);

    double x0 = 1.5;
    vd x = {0, 1, 2, 3, 4};
    vd y = {0, 1.8415, 2.9093, 3.1411, 3.2432};

    calculateSplineCoefficients(n, x, y, a, b, c, d);
    auto res = evaluate(n, x0, x, a, c, b, d);
    std::cout << "Function result in x0: " << res << std::endl;
}