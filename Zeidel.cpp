#include "src/Matrix.h"
 #define DEBUG;

using vd = std::vector<double>;

void simple_integrations(const matrix<double>& matrixA, const vd& b, const double& eps) {
    int n = matrixA.getCols();
    matrix<double> alpha(n, n);
    vd beta(n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i != j) {
                alpha[i][j] = (-1) * matrixA[i][j] / matrixA[i][i];
            }
        }
        beta[i] = b[i] / matrixA[i][i];
    }
#ifdef DEBUG
    std::cout << "Matrix alpha:\n" << alpha << "\n\n";
    std::cout << "Vector beta:\n";
    for (int i = 0; i < n; i++) {
        std::cout << beta[i] << '\n';
    }
#endif

    double norma = 0;
    for (int i = 0; i < n; i++) {
        double summ = 0;
        for (int j = 0; j < n; j++) {
            summ += fabs(alpha[i][j]);
        }
        if (norma < summ) {
            norma = summ;
        }
    }
#ifdef DEBUG
    std::cout << "Norma: " << norma << "\n\n";
#endif
    if (norma >= 1) {
        std::cout << "Еhe method of simple iterations "
                     "does not converge for this matrix.\n";
        return;
    }
    vd X_prev(beta);
    vd X_cur(beta);
    int k = 0;

    while (true) {
        k++;
        vd tmp = alpha * X_prev;
        for (int i = 0; i < n; i++) {
            X_cur[i] += tmp[i];
        }
        double norma_X = 0;
        for (int i = 0; i < n; i++) {
            if (double diff; (diff = fabs(X_cur[i] - X_prev[i])) > norma_X) {
                norma_X = diff;
            }
        }

#ifdef DEBUG
        std::cout << "X_prev    X_cur on iteration " << k << '\n';
        for (int i = 0; i < n; i++) {
            std::cout << X_prev[i] << "     " << X_cur[i] << '\n';
        }
        std::cout << "norma_X: " << norma_X * norma / (1 - norma) << '\n';
#endif

        if (fabs(norma_X * norma / (1 - norma)) <= eps) {
            break;
        }
        X_prev = X_cur;
        X_cur = beta;
    }
    std::cout << "Answer:\n";
    for (int i = 0; i < n; i++) {
        std::cout << std::fixed << std::setprecision(2) << X_cur[i] << '\n';
    }
    std::cout << "Number of iterations: " << k << '\n';
}

void Zeidel(const matrix<double>& matrixA, const vd& b, const double& eps) {
    int n = matrixA.getCols();
    matrix<double> alpha(n, n);
    vd beta(n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i != j) {
                alpha[i][j] = (-1) * matrixA[i][j] / matrixA[i][i];
            }
        }
        beta[i] = b[i] / matrixA[i][i];
    }
#ifdef DEBUG
    std::cout << "Matrix alpha:\n" << alpha << "\n\n";
    std::cout << "Vector beta:\n";
    for (int i = 0; i < n; i++) {
        std::cout << beta[i] << '\n';
    }
#endif

    double norma = 0;
    for (int i = 0; i < n; i++) {
        double summ = 0;
        for (int j = 0; j < n; j++) {
            summ += fabs(alpha[i][j]);
        }
        if (norma < summ) {
            norma = summ;
        }
    }
#ifdef DEBUG
    std::cout << "Norma: " << norma << "\n\n";
#endif
    if (norma >= 1) {
        std::cout << "Еhe method of simple iterations "
                     "does not converge for this matrix.\n";
        return;
    }
    vd X_prev(beta);
    vd X_cur(beta);
    int k = 0;

    while (true) {
        k++;
        for (int i = 0; i < n; i++) {
            X_cur[i] = beta[i];
            for (int j = 0; j < n; j++) {
                X_cur[i] += alpha[i][j] * X_cur[j];
            }
        }

        double norma_X = 0;
        for (int i = 0; i < n; i++) {
            if (double diff; (diff = fabs(X_cur[i] - X_prev[i])) > norma_X) {
                norma_X = diff;
            }
        }
        double norma_C = 0;
        for (int i = 0; i < n; i++) {
            double summ = 0;
            for (int j = i; j < n; j++) {
                summ += fabs(alpha[i][j]);
            }
            if (norma_C < summ) {
                norma_C = summ;
            }
        }

#ifdef DEBUG
        std::cout << "X_prev    X_cur on iteration " << k << '\n';
        for (int i = 0; i < n; i++) {
            std::cout << X_prev[i] << "     " << X_cur[i] << '\n';
        }
        std::cout << "norma_X: " << norma_X << '\n';
#endif

        if (fabs(norma_X) <= eps) {
            break;
        }
        X_prev = X_cur;
    }
    std::cout << "Answer:\n";
    for (int i = 0; i < n; i++) {
        std::cout << std::fixed << std::setprecision(2) << X_cur[i] << '\n';
    }
    std::cout << "Number of iterations: " << k << '\n';
}

int main() {
    constexpr double eps = 0.01;
    // const matrix<double> matrixA = {
    //     {10, 1, 1},
    //     {2, 10, 1},
    //     {2, 2, 10}
    // };
    // const vd b = {12, 13, 14};

    const matrix<double> matrixA = {
        {21, -6, -9, -4},
        {-6, 20, -4, 2},
        {-2, -7, -20, 3},
        {4, 9, 6, 24}
    };
    const vd b = {127, -144, 236, -5};

    simple_integrations(matrixA, b, eps);

    Zeidel(matrixA, b, eps);
}
