#include "src/Matrix.h"
#include <cmath>

//#define DEBUG

matrix<double> setE(size_t n);

int sign(double numb) {
    if (fabs(numb) < 0) {
        return 0;
    } else if (numb < 0) {
        return -1;
    }
    return 1;
}

matrix<double> setE(size_t n) {
    matrix<double> tmp(n, n);
    for (int i = 0; i < n; i++) {
        tmp[i][i] = 1;
    }
    return tmp;
}

void QR(const matrix<double>& matrixA, matrix<double>& Q, matrix<double>& R) {
    size_t n = matrixA.getCols();

    matrix<double> H(n, n);
    Q = setE(n);
    R = matrixA;

    for (int i = 0; i < n - 1; i++) {
        matrix<double> v(n, 1);
        for (int j = i; j < n; j++) {
            v[i][0] += R[j][i] * R[j][i];
            if (i != j) {
                v[j][0] = R[j][i];
            }
        }
        v[i][0] = R[i][i] + sign(R[i][i]) * sqrt(v[i][0]);
        matrix<double> v_t = v.transpose();
        double dev = 2 / (v_t * v)[0][0];
        H = setE(n);
        matrix<double> tmp = v * v_t;
        H -= (v * v_t) * dev;

        R = H * R;
        Q *= H;
#ifdef DEBUG
        std::cout << "Number of iteration: " << i << "\n";
        std::cout << "Vector V:\n";
        for (int k = 0; k < n; k++) {
            std::cout << v[k][0] << '\n';
        }
        std::cout << "Matrix H:\n" << H << "\n\nMatrix Q\n" << Q;
        std::cout << "\n\nMatrix R:\n" << R << "\n\nMatrix Q * R\n" << Q * R << "\n\n";
#endif
    }
#if defined (INFO) || defined(DEBUG)
    std::cout << "MatrixA:\n" << matrixA << "\n\nMatrix Q\n" << Q;
    std::cout << "\n\nMatrix R:\n" << R << "\n\nMatrix Q * R\n" << Q * R;
#endif
}

std::vector<std::complex<long double>> get_lambdas(const matrix<double> &matrixA, double eps, int max_iterations = 1000) {
    matrix<double> A(matrixA);
    if (A.getRows() != A.getCols()) {
        throw std::invalid_argument("Matrix must be square for QR algorithm");
    }

    int n = A.getRows();
    std::vector<std::complex<long double>> lambdas(n, 0.0);

    int iteration = 0;

    while (iteration < max_iterations) {
        iteration++;
        matrix<double> Q = setE(n);
        matrix<double> R = setE(n);

        QR(A, Q, R);
        A = R * Q;

        bool isLowZero = true;
        for (int i = 1; i < n; i++) {
            for (int j = 0; j < i; j++) {
                if (std::abs(A[i][j]) > eps) {
                    isLowZero = false;
                    break;
                }
            }
            if (!isLowZero) {
                break;
            }
        }

        if (isLowZero) {
            break;
        }
    }

    std::cout << "QR algorithm completed in " << iteration << " iterations" << std::endl;
    if (iteration >= max_iterations) {
        std::cout << "Warning: Maximum iterations reached" << std::endl;
    }

    for (int i = 0; i < n; i++) {
        if (i < n - 1 && std::abs(A[i + 1][i]) > eps) {
            long double a = A[i][i];
            long double b = A[i][i + 1];
            long double c = A[i + 1][i];
            long double d = A[i + 1][i + 1];

            long double trace = a + d;
            long double determinant = a * d - b * c;
            long double discriminant = trace * trace - 4 * determinant;

            if (discriminant < 0) {
                long double real_part = trace / 2;
                long double imag_part = std::sqrt(-discriminant) / 2;
                lambdas[i] = {real_part, imag_part};
                lambdas[i + 1] = {real_part, -imag_part};
                i++;
            } else {
                lambdas[i] = A[i][i];
            }
        } else {
            lambdas[i] = A[i][i];
        }
    }
    return lambdas;
}

int main() {
    constexpr double eps = 0.01;
//     const matrix<double> matrixA = {
//         {1, 3, 1},
//         {1, 1, 4},
//         {4, 3, 1}
//     };

     const matrix<double> matrixA = {
         {-9, 2, 2},
         {-2, 0, 7},
         {8, 2, 0}
     };


    matrix<double> Q, R;
    QR(matrixA, Q, R);

    std::vector<std::complex<long double>> lambdas = get_lambdas(matrixA, 0.001, 1000);
    std::cout << "Lambdas:" << "\n";
    for (auto lambda : lambdas) {
        std::cout << std::fixed <<std::setprecision(2) << lambda << std::endl;
    }
}

