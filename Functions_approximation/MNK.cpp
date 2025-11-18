#include <iostream>
#include "../LU_Decomposition.cpp"

using vd = std::vector<double>;


vd sub_first_approximation(const vd& x, const vd& y) {
    matrix<double> A(2, 2);
    vd b(2);
    A[0][0] = x.size();
    for (int i = 0; i < x.size(); ++i) {
        A[0][1] += x[i];
        A[1][1] += x[i] * x[i];
        b[0] += y[i];
        b[1] += x[i] * y[i];
    }
    A[1][0] = A[0][1];

    matrix<double> L(2, 2);
    matrix<double> U(2, 2);
    vd X(2);
    splitLU(2, A, L, U);
    solve(2, L, U, b, X);
    return X;
}

double first_approximation(const double x0, const vd& X) {
    return X[0] + X[1] * x0;
}

vd sub_second_approximation(const vd& x, const vd& y) {
    matrix<double> A(3, 3);
    vd b(3);
    A[0][0] = x.size();
    for (int i = 0; i < x.size(); ++i) {
        A[0][1] += x[i];
        A[1][1] += x[i] * x[i];
        A[1][2] += x[i] * x[i] * x[i];
        A[2][2] += x[i] * x[i] * x[i] * x[i];

        b[0] += y[i];
        b[1] += x[i] * y[i];
        b[2] += x[i] * x[i] * y[i];
    }
    A[1][0] = A[0][1];
    A[2][0] = A[1][1];
    A[0][2] = A[1][1];
    A[2][1] = A[1][2];

    matrix<double> L(3, 3);
    matrix<double> U(3, 3);
    vd X(3);
    splitLU(3, A, L, U);
    solve(3, L, U, b, X);
    return X;
}

double second_approximation(const double x0, const vd& X) {
    return X[0] + X[1] * x0 + X[2] * x0 * x0;
}

double first_sum_of_squares_errors(const vd& x, const vd& y) {
    double F = 0;
    const vd A = sub_first_approximation(x, y);

    for (int i = 0; i < x.size(); ++i) {
        const double tmp = first_approximation(x[i], A) - y[i];
        F += tmp * tmp;
    }
    return F;
}

double second_sum_of_squares_errors(const vd& x, const vd& y) {
    double F = 0;
    const vd A = sub_second_approximation(x, y);

    for (int i = 0; i < x.size(); ++i) {
        const double tmp = second_approximation(x[i], A) - y[i];
        F += tmp * tmp;
    }
    return F;
}

int main() {
    // vd x = {0, 1.7, 3.4, 5.1, 6.8, 8.5}; // example
    // vd y = {0, 1.3038, 1.8439, 2.2583, 2.6077, 2.9155};

    vd x = {-0.7, -0.4, -0.1, 0.2, 0.5, 0.8};
    vd y = {-0.7754, -0.41152, -0.10017, 0.20136, 0.5236, 0.9273};

    std::cout << std::fixed << std::setprecision(5);
    const vd A = sub_first_approximation(x, y);
    std::cout << "First approximation:\nx   F1" << std::endl;
    for (const double i : x) {
        std:: cout << i << "  " << first_approximation(i, A) << '\n';
    }
    for (int i = 0; i < A.size(); ++i) {
        std::cout << 'a' << i << " = " << A[i] << '\n';
    }
    std::cout << "First sum of squares errors: " << first_sum_of_squares_errors(x, y) << '\n';

    const vd B = sub_second_approximation(x, y);
    std::cout << "First approximation:\nx   F2" << std::endl;
    for (const double i : x) {
        std:: cout << i << "  " << second_approximation(i, B) << '\n';
    }
    for (int i = 0; i < B.size(); ++i) {
        std::cout << 'a' << i << " = " << B[i] << '\n';
    }
    std::cout << "Second sum of squares errors: " << second_sum_of_squares_errors(x, y) << '\n';
}