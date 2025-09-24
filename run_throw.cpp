#include <bits/stdc++.h>


using vvd = std::vector<std::vector<double>>;
using vd = std::vector<double>;


int main() {
    int n = 5;
//    vvd matrixA = {{8, -2, 0, 0},
//                   {-1, 6, -2, 0},
//                   {0, 2, 10, -4},
//                   {0, 0, -1, 6}};
//    vd d = {6, 3, 8, 5};
    vvd matrixA = {{13, -5, 0, 0, 0},
                   {-4, 9, -5, 0, 0, 0},
                   {0, -1, 12, -6, 0},
                   {0, 0, 6, 20, -5},
                   {0, 0, 0, 4, 5}};
    vd d = {-66, -47, -43, -74, 14};
    vd P(n);
    vd Q(n);
    vd X(n);
    P[0] = (-1) * matrixA[0][1]/matrixA[0][0];
    Q[0] = d[0]/matrixA[0][0];

    for (int i = 1; i < n; ++i) {
        double denom = matrixA[i][i] + matrixA[i][i - 1] * P[i - 1];
        P[i] = (-1) * matrixA[i][i + 1] / denom;
        Q[i] = (d[i] - matrixA[i][i - 1] * Q[i - 1]) / denom;
    }

    X[n - 1] = (d[n - 1] - matrixA[n - 1][n - 2] * Q[n - 2])
            / (matrixA[n - 1][n - 1] + matrixA[n - 1][n - 2] * P[n - 2]);

    X[n - 1] = Q[n - 1];
    for (int i = n - 2; i >= 0; i--) {
        X[i] = P[i] * X[i + 1] + Q[i];
    }
    for (int i = 0; i < n; ++i) {
        std::cout << X[i] << '\n';
    }
}
