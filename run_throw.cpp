#include <bits/stdc++.h>


using vvd = std::vector<std::vector<double>>;
using vd = std::vector<double>;


int main() {
    int n = 5;
    vvd matrixA = {{-11, -8, 0, 0, 0},
                   {9, -17, 1, 0, 0, 0},
                   {0, -4, 20, 9, 0},
                   {0, 0, -4, 14, 3},
                   {0, 0, 0, -6, 14}};
    vd d = {99, -75, 66, 54, 8};
    vd P(n);
    vd Q(n);
    vd X(n);
    P[0] = -matrixA[0][1]/matrixA[0][0];
    Q[0] = d[0]/matrixA[0][0];

    for (int i = 1; i < n; ++i) {
        double denom = matrixA[i][i] + matrixA[i][i - 1] * P[i - 1];
        P[i] = -matrixA[i][i + 1] / denom;
        Q[i] = (d[i] - matrixA[i][i - 1] * Q[i - 1]) / denom;
    }

//    X[n - 1] = (d[n - 1] - matrixA[n - 2][n - 1] * Q[n - 2])
//            / (matrixA[n - 1][n - 1] + matrixA[n - 2][n - 1] * P[n - 1]);

    X[n - 1] = Q[n - 1];
    for (int i = n - 2; i >= 0; --i) {
        X[i] = P[i] * X[i + 1] + Q[i];
    }
    for (int i = 0; i < n; ++i) {
        std::cout << X[i] << '\n';
    }
}
