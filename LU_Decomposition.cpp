#include <cmath>
#include <iostream>
#include <math.h>
#include <vector>
#include "src/Matrix.h"

using vvd = matrix<double>;
using vd = std::vector<double>;


//bool check_max_diag(const int n, const vvd &matrix)  {// TODO: как поравнять и в строках и в столбцах
//    for (int i = 0; i < n; ++i) {
//        int summ = 0;
//
//        for (int j = 0; j < n; ++j) {
//            summ += matrix[i][j];
//        }
//        if (matrix[i][i] <= summ) {
//            return false;
//        }
//    }
//    return true;
//}

enum {
    SUCCESS,
    INVALID_MATRIX,
    ZERO_DET,
    NO_DIAG_DOMINATION,
    ERROR
};

void ValidateCode(const int code) {
    switch (code) {
        case SUCCESS:
            printf("SUCCESS\n");
        return;
        case ERROR:
            printf("ERROR\n");
        return;
        case INVALID_MATRIX:
            printf("First element cant be zero.\n");
        return;
        case NO_DIAG_DOMINATION:
            printf("There is no diagonal domination.\n");
        return;
        case ZERO_DET:
            printf("There is no single solution.\n");
        return;
        default:
            printf("why are u here?\n");
        return;
    }
}

// double matrixDet(const int n, const vvd &matrix) { //в тупую
//     switch (n) {
//         case 1:
//             return matrix[0][0];
//         case 2:
//             return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
//         case 3: { // метод треугольников
//             const double summ = matrix[0][0] * matrix[1][1] * matrix[2][2] + matrix[0][1] * matrix[1][2] * matrix[2][0]
//             + matrix[1][0] * matrix[2][1] * matrix[0][2];
//             const double diff = matrix[0][2] * matrix[1][1] * matrix[2][0] + matrix[1][0] * matrix[0][1] * matrix[2][2]
//             + matrix[1][2] * matrix[2][1] * matrix[0][0];
//             return summ - diff;
//         }
//         default:
//             if (n < 1) {
//                 std::cout << "Something went wrong in matrixDet\n";
//                 throw ERROR;
//             }
//             break;
//     }
//     double det = 0;
//
//     for (int i = 0; i < n; ++i) {
//         vvd temp(matrix);
//         temp.pop_back();
//         for (int j = 0; j < n - 1; ++j) {
//             temp[j].erase(temp[j].begin() + i);
//         }
//         det += matrix[n - 1][i] * ((i + n + 1) % 2 == 0 ? 1 : -1) * matrixDet(n - 1, temp);
//     }
//     return det;
// }

void splitLU(const int n, const vvd &matrix, vvd &L, vvd &U) {
    for (int k = 0; k < n; ++k) {
        U[0][k] = matrix[0][k];
    }
    L[0][0] = 1;
    for (int k = 1; k < n; ++k) {
        L[k][0] = matrix[k][0] / U[0][0];
    }

    //первые 3 строчки и столбца по отдельности
    // for (int k = 2; k < n; ++k) {
    //     U[2][k] = matrix[2][k] - L[2][1] * U[1][k];
    // }
    // for (int k = 3; k < n; ++k) {
    //     L[k][2] = (matrix[k][2] - L[k][1] * U[1][2]) / U[2][2];
    // }
    // for (int k = 3; k < n; ++k) {
    //     U[3i][k] = matrix[3i][k] - L[3i][1] * U[1][k] - L[3i][2] * U[2][k];
    // }
    // for (int k = 4; k < n; ++k) {
    //     L[k][3i] = (matrix[k][3i] - L[k][1] * U[1][3i] - L[k][2] * U[2][3i]) / U[3][3];
    // }

    //все елементы кроме 1 строчки и столбца в общем виде
    for (int i = 1; i < n; ++i){
        for (int k = i; k < n; ++k) {
            U[i][k] = matrix[i][k];
            for (int j = 0; j < i; ++j) {
                U[i][k] -= L[i][j] * U[j][k];
                if (k != i) {
                    L[k][i] -= L[k][j] * U[j][i];
                }
            }
            if (k != i) {
                L[k][i] += matrix[k][i];
                L[k][i] /= U[i][i];
            } else {
                L[k][i] = 1;
            }
        }
    }
}

void solve(const int n, const vvd &L, const vvd &U, const vd &b, vd &x) {
    //Lz = b
    vd z(n);
    z[0] = b[0];
    for (int i = 1; i < n; ++i) {
        z[i] = b[i];
        for (int j = 0; j < i; ++j) {
            z[i] -= L[i][j] * z[j];
        }
    }

    //Ux = z
    x[n - 1] = z[n - 1] / U[n - 1][n - 1];
    for (int i = n - 2; i >= 0; --i) {
        x[i] = z[i];
        for (int j = i + 1; j < n; ++j) {
            x[i] -= U[i][j] * x[j];
        }
        x[i] /= U[i][i];
    }
}

int inverseMatrix(const int n, const vvd &matrixA, vvd &Ansv) {
    vvd L(n, n), U(n, n);
    splitLU(n, matrixA, L, U);

    for (int i = 0; i < n; ++i) {
        vd b(n);
        vd x(n);
        b[i] = 1;
        solve(n, L, U, b, x);

        for (int j = 0; j < n; ++j) {
            Ansv[j][i] =x[j];
        }
    }
    return SUCCESS;
}

double DetLU(const int n, const vvd &matrixA) {
    vvd L(n, n);
    vvd U(n, n);
    splitLU(n, matrixA, L, U);
    double det = 1;
    for (int i = 0; i < n; ++i) {
        det *= U[i][i];
    }
    return det;
}

int main() {
    int n;
    std::cout << "Input number of variables:\n";
    std::cin >> n;

    // vvd matrixA = {{10, 1, 1},
    //                {2, 10, 1},
    //                {2, 2, 10}};
    // vd b = {12, 13, 14};

    // vvd matrixA = {{-4, -9, 4, 3},
    //                {2, 7, 9, 8},
    //                {4, -4, 0, -2},
    //                {-8, 5, 2, 9}};
    // vd b = {-52, 76, 26, -73};

    // vvd matrixA = {{-11, -8, 0, 0, 0},
    //                {9, -17, 1, 0, 0, 0},
    //                {0, -4, 20, 9, 0},
    //                {0, 0, -4, 14, 3},
    //                {0, 0, 0, -6, 14}};
    // vd b = {99, -75, 66, 54, 8};

    vvd matrixA = {{13, -5, 0, 0, 0},
                   {-4, 9, -5, 0, 0},
                   {0, -1, 12, -6, 0},
                   {0, 0, 6, 20, -5},
                   {0, 0, 0, 4, 5}};
    const vd b = {-66, -47, -43, -74, 14};

    if (fabs(matrixA[0][0]) <= 0.0000000001) {
        return INVALID_MATRIX;
    }
    vd X(n);
    vvd L(n, n);
    vvd U(n, n);
    splitLU(n, matrixA, L, U);

    std::cout << "Matrix L:\n" << L << "\n\nMatrix U:\n" << U;
    std::cout << "\n\nL * U:\n" << L * U << "\n";

    solve(n, L, U, b, X);

    std::cout << "\nAnswers:\n";
    for (int i = 0; i < n; ++i) {
        std::cout << X[i] << "\n";
    }

    std::cout << "\nInvers matrix: \n";
    vvd Inv(n, n);
    inverseMatrix(n, matrixA, Inv);
    std::cout << Inv << '\n';

    double determenant = DetLU(n, matrixA);
    std::cout << "\nDetermenant: " << determenant << '\n';
}