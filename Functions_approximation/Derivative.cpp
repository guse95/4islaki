#include <iostream>
#include <vector>

using vd = std::vector<double>;

double first_derivative(const double x0, const vd& x, const vd& y) {
    int segment = 0;
    for (int i = 0; i <= x.size(); i++) {
        if (x0 <= x[i]) {
            segment = i - 1;
            break;
        }
    }

    double res = (y[segment + 1] - y[segment]) /
        (x[segment + 1] - x[segment]);

    double tmp = (y[segment + 2] - y[segment + 1]) /
        (x[segment + 2] - x[segment + 1]);

    tmp -= res;
    tmp /= x[segment + 2] - x[segment];
    tmp *= (2 * x0 - x[segment] - x[segment + 1]);
    res += tmp;
    return res;
}

double second_derivative(const double x0, const vd& x, const vd& y) {
    int segment = 0;
    for (int i = 0; i <= x.size(); i++) {
        if (x0 <= x[i]) {
            segment = i - 1;
            break;
        }
    }

    double tmp = (y[segment + 1] - y[segment]) /
        (x[segment + 1] - x[segment]);

    double res = (y[segment + 2] - y[segment + 1]) /
        (x[segment + 2] - x[segment + 1]);

    res -= tmp;
    res /= (x[segment + 2] - x[segment]) / 2;
    return res;
}

int main() {
    double x0 = 0.2;
    vd x = {0.0, 0.1, 0.2, 0.3, 0.4};
    vd y = {1.0, 1.1052, 1.2214, 1.3499, 1.4918};

    std::cout << "First derivative: " << first_derivative(x0, x, y) << std::endl;
    std::cout << "Second derivative: " << second_derivative(x0, x, y) << std::endl;
}