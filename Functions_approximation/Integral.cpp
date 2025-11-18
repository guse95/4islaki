#include <iomanip>
#include <iostream>
#include <vector>
#include <cmath>

using vd = std::vector<double>;


double f(const double x) {
    // return x / ((3 * x + 4) * (3 * x + 4));
    return 1 / (x * x + 4);
}

double F(const double x) {
    return std::atan(x / 2) / 2;
    // return std::log(fabs(3 * x + 4)) / 9 + 4 / (9 * (3 * x + 4));
}

double rectangles(const double l, const double r, const double h) {
    double res = 0;
    double x_i = l;

    while (x_i + h < r) {
        res += f((2 * x_i + h) / 2);
        x_i += h;
    }

    res += f((x_i + r) / 2);
    res *= h;
    return res;
}

double trapezoids(const double l, const double r, const double h) {
    double res = f(l) / 2;
    double x_i = l + h;

    while (x_i < r) {
        res += f(x_i);
        x_i += h;
    }

    res += f(r) / 2;
    res *= h;
    return res;
}

double Simpson(const double l, const double r, const double h) {
    double res = f(l);
    double x_i = l + h;
    int coef = 4;

    while (x_i < r) {
        res += f(x_i) * coef;
        x_i += h;
        coef = (coef == 2) ? 4 : 2;
    }

    res += f(r);
    res *= h;
    res /= 3;
    return res;
}

int main() {
    // constexpr double x0 = -1;
    // constexpr double x1 = 1;
    // constexpr double h1 = 0.5;
    // constexpr double h2 = 0.25;

    constexpr double x0 = -2;
    constexpr double x1 = 2;
    constexpr double h1 = 1;
    constexpr double h2 = 0.5;


    double F_rec_1 = rectangles(x0, x1, h1);
    double F_trap_1 = trapezoids(x0, x1, h1);
    double F_sim_1 = Simpson(x0, x1, h1);
    std::cout << std::fixed << std::setprecision(5);
    std::cout << "Padding: " << h1 << std::endl;
    std:: cout << "Rectangles method: " << F_rec_1 << std::endl;
    std:: cout << "Trapezoids method: " << F_trap_1 << std::endl;
    std:: cout << "Simpson`s method: " << F_sim_1 << std::endl;

    double F_rec_2 = rectangles(x0, x1, h2);
    double F_trap_2 = trapezoids(x0, x1, h2);
    double F_sim_2 = Simpson(x0, x1, h2);
    std::cout << "Padding: " << h2 << std::endl;
    std:: cout << "Rectangles method: " << F_rec_2 << std::endl;
    std:: cout << "Trapezoids method: " << F_trap_2 << std::endl;
    std:: cout << "Simpson`s method: " << F_sim_2 << std::endl;

    std::cout << "Runge-Romberg-Richardson integral:" << std::endl;
    double k = h1 / h2;
    double F_rec_specified = F_rec_1 + (F_rec_1 - F_rec_2) /
        (k * k - 1);
    double F_trap_specified = F_trap_1 + (F_trap_1 - F_trap_2) /
        (k * k - 1);
    double F_sim_specified = F_sim_1 + (F_sim_1 - F_sim_2) /
        (k * k * k * k - 1);

    double F_accurate = F(x1);
    F_accurate -= F(x0);
    std::cout << "Accurate integral: " << F_accurate << std::endl;
    std:: cout << "Rectangles method: " << F_rec_specified << "; absolute inaccuracy = " <<
        fabs(F_accurate - F_rec_specified) << std::endl;
    std:: cout << "Trapezoids method: " << F_trap_specified << "; absolute inaccuracy = " <<
        fabs(F_accurate - F_trap_specified) << std::endl;
    std:: cout << "Simpson`s method: " << F_sim_specified << "; absolute inaccuracy = " <<
        fabs(F_accurate - F_sim_specified) << std::endl;
}