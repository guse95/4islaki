#include <cmath>
#include <iostream>
#include <math.h>
#include <vector>

using vvd = std::vector<std::vector<double>>;
using vd = std::vector<double>;

///Не нужно в 1 задаче
bool check_max_diag(const int n, const vvd &matrix)  {// TODO: wtf как поравнять и в строках и в столбцах
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i != j && matrix[i][i] <= matrix[i][j]) {
                return false;
            }
        }
    }
    return true;
}

enum {
    SUCCESS,
    INVALID_MATRIX,
    ERROR
};

double matrixDet(const int n, const vvd &matrix) { //в тупую
    switch (n) {
        case 1:
            return matrix[0][0];
        case 2:
            return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
        case 3: { // метод треугольников
            const double summ = matrix[0][0] * matrix[1][1] * matrix[2][2] + matrix[0][1] * matrix[1][2] * matrix[2][0]
            + matrix[1][0] * matrix[2][1] * matrix[0][2];
            const double diff = matrix[0][2] * matrix[1][1] * matrix[2][0] + matrix[1][0] * matrix[0][1] * matrix[2][2]
            + matrix[1][2] * matrix[2][1] * matrix[0][0];
            return summ - diff;
        }
        default:
            if (n < 1) {
                std::cout << "Something went wrong in matrixDet\n";
                throw ERROR;
            }
            break;
    }
    double det = 0;

    for (int i = 0; i < n; ++i) {
        vvd temp(matrix);
        temp.pop_back();
        for (int j = 0; j < n - 1; ++j) {
            temp[j].erase(temp[j].begin() + i);
        }
        det += matrix[n - 1][i] * ((i + n + 1) % 2 == 0 ? 1 : -1) * matrixDet(n - 1, temp);
    }
    return det;
}

void splitLU(const int n, const vvd &matrix, vvd &L, vvd &U) {
    for (int k = 0; k < n; ++k) {
        U[1][k] = matrix[1][k];
    }
    for (int k = 1; k < n; ++k) {
        L[k][1] = matrix[k][1] / U[1][k];
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
            for (int j = 1; j < i; ++j) {
                U[i][k] -= L[i][j] * U[j][k];
                if (k != i) {
                    L[k][i] = L[k][j] * U[j][i];
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

int solve(const int n, const vvd &matrixA, const vd &b, vd &X) {
    if (matrixA[0][0] == 0 || matrixDet(n, matrixA) == 0) {
        return INVALID_MATRIX;
    }

    vvd L(n), U(n);
    splitLU(n, matrixA, L, U);

}

int main() {
    int n;
    double elem;
    std::cout << "Input number of variables:\n";
    std::cin >> n;
    std::cout << "Input matrix:\n";
    vvd matrixA(n, vd(n));
    vd b(n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            std::cin >> elem;
            matrixA[i][j] = elem;
        }
        std::cin >> elem;
        b[i] = elem;
    }

//    check_max_diag(n, matrixA);
    //TODO: функция по разложению на матрицы LU
}
