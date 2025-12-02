#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <functional>

using namespace std;

struct Solution {
    vector<double> x;
    vector<double> y;
};

struct SystemSolution {
    vector<double> x;
    vector<double> y;
    vector<double> z; // y' = z
};

double exact_solution(double x) {
    return 3*x + exp(-x*x);
}

auto f_ode = [](double x, double y, double z) {
    if (abs(2*x + 1) < 1e-10) {
        return 0.0; // Handle singularity
    }
    return (4*y - 4*x*z) / (2*x + 1);
};

SystemSolution runge_kutta_system(double x0, double y0, double z0, double h, int n,
                                  function<double(double, double, double)> f,
                                  function<double(double, double, double)> g) {
    SystemSolution sol;
    sol.x.resize(n + 1);
    sol.y.resize(n + 1);
    sol.z.resize(n + 1);

    sol.x[0] = x0;
    sol.y[0] = y0;
    sol.z[0] = z0;

    for (int k = 0; k < n; k++) {
        double xk = sol.x[k];
        double yk = sol.y[k];
        double zk = sol.z[k];

        double K1 = h * f(xk, yk, zk);
        double L1 = h * g(xk, yk, zk);

        double K2 = h * f(xk + h/2, yk + K1/2, zk + L1/2);
        double L2 = h * g(xk + h/2, yk + K1/2, zk + L1/2);

        double K3 = h * f(xk + h/2, yk + K2/2, zk + L2/2);
        double L3 = h * g(xk + h/2, yk + K2/2, zk + L2/2);

        double K4 = h * f(xk + h, yk + K3, zk + L3);
        double L4 = h * g(xk + h, yk + K3, zk + L3);

        sol.x[k + 1] = xk + h;
        sol.y[k + 1] = yk + (K1 + 2*K2 + 2*K3 + K4) / 6;
        sol.z[k + 1] = zk + (L1 + 2*L2 + 2*L3 + L4) / 6;
    }

    return sol;
}

Solution shooting_method_mixed(double a, double b,
                              function<double(double, double, double)> boundary_condition_a,
                              double boundary_value_a,
                              function<double(double, double, double)> boundary_condition_b,
                              double boundary_value_b,
                              function<double(double, double, double)> f_ode,
                              double eta0, double eta1, double h, double epsilon) {
    auto f = [](double x, double y, double z) { return z; }; // y' = z
    auto g = f_ode; // z' = f_ode(x, y, z)

    vector<double> eta_vals = {eta0, eta1};
    vector<double> phi_vals;

    for (int i = 0; i < 2; i++) {
        SystemSolution sol = runge_kutta_system(a, 0, eta_vals[i], h,
                                               int((b - a) / h), f, g);

        double y_a = sol.y[0];
        double z_a = sol.z[0];
        double phi = boundary_condition_a(y_a, z_a, a) - boundary_value_a;
        phi_vals.push_back(phi);
    }

    int iteration = 0;
    while (abs(phi_vals.back()) > epsilon && iteration < 100) {
        double eta_next = eta_vals[eta_vals.size() - 1] -
                         (eta_vals[eta_vals.size() - 1] - eta_vals[eta_vals.size() - 2]) /
                         (phi_vals[phi_vals.size() - 1] - phi_vals[phi_vals.size() - 2]) *
                         phi_vals[phi_vals.size() - 1];

        SystemSolution sol = runge_kutta_system(a, 0, eta_next, h,
                                               int((b - a) / h), f, g);

        double y_a = sol.y[0];
        double z_a = sol.z[0];
        double phi_next = boundary_condition_a(y_a, z_a, a) - boundary_value_a;

        eta_vals.push_back(eta_next);
        phi_vals.push_back(phi_next);
        iteration++;
    }

    SystemSolution sol_final = runge_kutta_system(a, 0, eta_vals.back(), h,
                                                 int((b - a) / h), f, g);

    double z_at_0 = 0;
    double y_at_0 = 0;
    for (size_t i = 0; i < sol_final.x.size(); i++) {
        if (abs(sol_final.x[i]) < 1e-10) {
            z_at_0 = sol_final.z[i];
            y_at_0 = sol_final.y[i];
            break;
        }
    }

    double scale_factor = 1.0 / z_at_0;
    double constant_C = -y_at_0 * scale_factor;

    Solution result;
    result.x = sol_final.x;
    result.y.resize(sol_final.y.size());

    for (size_t i = 0; i < sol_final.y.size(); i++) {
        result.y[i] = constant_C + scale_factor * sol_final.y[i];
    }

    return result;
}

Solution finite_difference_mixed(double a, double b,
                                function<double(double)> p,
                                function<double(double)> q,
                                function<double(double)> f_func,
                                double h,
                                function<double(double, double)> bc_left,
                                function<double(double, double)> bc_right) {
    int n = int((b - a) / h);
    vector<double> x(n + 1);
    vector<double> y(n + 1);

    for (int i = 0; i <= n; i++) {
        x[i] = a + i * h;
    }

    vector<double> A(n + 1), B(n + 1), C(n + 1), D(n + 1);
    A[0] = 2 - 1/h;
    B[0] = 1/h;
    C[0] = 0;
    D[0] = -9;

    for (int i = 1; i < n; i++) {
        double xi = x[i];
        double p_val = 4*xi/(2*xi + 1);
        double q_val = -4/(2*xi + 1);

        A[i] = 1/h/h - p_val/(2*h);
        B[i] = -2/h/h + q_val;
        C[i] = 1/h/h + p_val/(2*h);
        D[i] = 0;
    }

    A[n] = -1/h;
    B[n] = 1/h;
    C[n] = 0;
    D[n] = 1;

    vector<double> alpha(n + 1), beta(n + 1);
    alpha[0] = -B[0] / A[0];
    beta[0] = D[0] / A[0];

    for (int i = 1; i <= n; i++) {
        double denom = A[i] + C[i] * alpha[i-1];
        alpha[i] = -B[i] / denom;
        beta[i] = (D[i] - C[i] * beta[i-1]) / denom;
    }

    y[n] = beta[n];
    for (int i = n-1; i >= 0; i--) {
        y[i] = alpha[i] * y[i+1] + beta[i];
    }

    Solution sol;
    sol.x = x;
    sol.y = y;
    return sol;
}

double runge_romberg_error(const Solution& sol_h, const Solution& sol_2h, int p) {
    int n = sol_h.x.size() - 1;
    int mid_idx = n/2;
    double y_h = sol_h.y[mid_idx];

    double x_mid = sol_h.x[mid_idx];
    double y_2h = 0;
    for (size_t i = 0; i < sol_2h.x.size(); i++) {
        if (abs(sol_2h.x[i] - x_mid) < 1e-10) {
            y_2h = sol_2h.y[i];
            break;
        }
    }

    return abs(y_h - y_2h) / (pow(2, p) - 1);
}

void print_solution_table(const Solution& sol, const string& method_name) {
    cout << "\nTable. Solution by " << method_name << " method\n";
    cout << "=====================================================\n";
    cout << setw(8) << "x_k" << setw(15) << "y_k"
         << setw(15) << "y_exact" << setw(15) << "Error" << endl;
    cout << "=====================================================\n";

    for (size_t i = 0; i < sol.x.size(); i++) {
        double y_exact = exact_solution(sol.x[i]);
        double error = abs(y_exact - sol.y[i]);
        cout << fixed << setprecision(5) << setw(8) << sol.x[i]
             << setw(15) << sol.y[i]
             << setw(15) << y_exact
             << setw(15) << setprecision(6) << error << endl;
    }
    cout << "=====================================================\n";
}

void print_comparison_table(const Solution& sol1, const Solution& sol2,
                           const string& name1, const string& name2) {
    cout << "\nComparison of " << name1 << " and " << name2 << " methods\n";
    cout << "===============================================================\n";
    cout << setw(8) << "x_k" << setw(20) << name1
         << setw(20) << name2 << setw(15) << "Difference" << endl;
    cout << "===============================================================\n";

    for (size_t i = 0; i < sol1.x.size(); i++) {
        double diff = abs(sol1.y[i] - sol2.y[i]);
        cout << fixed << setprecision(5) << setw(8) << sol1.x[i]
             << setw(20) << sol1.y[i]
             << setw(20) << sol2.y[i]
             << setw(15) << setprecision(6) << diff << endl;
    }
    cout << "===============================================================\n";
}

int main() {
    cout << "==============================================\n";
    cout << "   BOUNDARY VALUE PROBLEM SOLUTION\n";
    cout << "   (2x+1)y'' + 4xy' - 4y = 0\n";
    cout << "   y'(-2) + 2y(-2) = -9\n";
    cout << "   y'(0) = 1\n";
    cout << "   Exact solution: y(x) = 3x + e^(-x^2)\n";
    cout << "   Interval: [-2, 0]\n";
    cout << "==============================================\n";

    double a = -2.0, b = 0.0;
    double h = 0.2;
    double epsilon = 1e-6;

    auto bc_left = [](double y, double z, double x) {
        return z + 2*y;
    };
    double bc_left_value = -9;

    auto bc_right_condition = [](double y, double z, double x) {
        return z; // y'
    };
    double bc_right_value = 1;

    cout << "\n1. SHOOTING METHOD\n";

    double eta0 = 1.0;
    double eta1 = 2.0;

    Solution sol_shooting = shooting_method_mixed(a, b, bc_left, bc_left_value,
                                                 bc_right_condition, bc_right_value,
                                                 f_ode, eta0, eta1, h, epsilon);

    print_solution_table(sol_shooting, "shooting");

    cout << "\n2. FINITE DIFFERENCE METHOD\n";

    auto p_func = [](double x) {
        return 4*x/(2*x + 1);
    };

    auto q_func = [](double x) {
        return -4/(2*x + 1);
    };

    auto f_func = [](double x) {
        return 0.0;
    };

    auto bc_left_fd = [h](double y0, double y1) {
        return (y1 - y0)/h + 2*y0; // y' + 2y
    };

    auto bc_right_fd = [h](double yn, double yn1) {
        return (yn - yn1)/h; // y'
    };

    Solution sol_fd = finite_difference_mixed(a, b, p_func, q_func, f_func, h,
                                            bc_left_fd, bc_right_fd);

    print_solution_table(sol_fd, "finite difference");

    cout << "\n3. ERROR ESTIMATION USING RUNGE-ROMBERG METHOD\n";

    Solution sol_shooting_h = shooting_method_mixed(a, b, bc_left, bc_left_value,
                                                   bc_right_condition, bc_right_value,
                                                   f_ode, eta0, eta1, h, epsilon);

    Solution sol_shooting_h2 = shooting_method_mixed(a, b, bc_left, bc_left_value,
                                                    bc_right_condition, bc_right_value,
                                                    f_ode, eta0, eta1, h/2, epsilon*2);

    int p = 4;
    double error_RR = runge_romberg_error(sol_shooting_h, sol_shooting_h2, p);

    cout << "Step h = " << h << endl;
    cout << "Step h/2 = " << h/2 << endl;
    cout << "Method order p = " << p << endl;
    cout << "Error estimate by Runge-Romberg: "
    << setprecision(6) << error_RR << endl;

    double x_mid = (a + b) / 2;
    double y_exact_mid = exact_solution(x_mid);
    double y_numeric_mid = sol_shooting_h.y[sol_shooting_h.x.size()/2];
    double actual_error = abs(y_exact_mid - y_numeric_mid);

    cout << "Actual error at x = " << x_mid << ": "
    << setprecision(6) << actual_error << endl;

    cout << "\n4. COMPARISON OF NUMERICAL METHODS\n";
    print_comparison_table(sol_shooting, sol_fd, "Shooting", "Finite Difference");

    cout << "\n5. SUMMARY OF RESULTS AT SELECTED POINTS\n";
    cout << "===================================================================================\n";
    cout << setw(8) << "x" << setw(15) << "Exact" << setw(15) << "Shooting"
         << setw(15) << "Finite Diff" << setw(15) << "Error Shoot"
         << setw(15) << "Error FD" << endl;
    cout << "===================================================================================\n";

    vector<double> test_points = {-2.0, -1.5, -1.0, -0.5, 0.0};
    for (double x : test_points) {
        double y_exact = exact_solution(x);

        double y_shoot = 0;
        for (size_t i = 0; i < sol_shooting.x.size(); i++) {
            if (abs(sol_shooting.x[i] - x) < 1e-10) {
                y_shoot = sol_shooting.y[i];
                break;
            }
        }

        double y_fd = 0;
        for (size_t i = 0; i < sol_fd.x.size(); i++) {
            if (abs(sol_fd.x[i] - x) < 1e-10) {
                y_fd = sol_fd.y[i];
                break;
            }
        }

        double error_shoot = abs(y_exact - y_shoot);
        double error_fd = abs(y_exact - y_fd);

        cout << fixed << setprecision(3) << setw(8) << x
        << setprecision(8) << setw(15) << y_exact
        << setw(15) << y_shoot
        << setw(15) << y_fd
        << setprecision(6) << setw(15) << error_shoot
        << setw(15) << error_fd << endl;
    }
    cout << "===================================================================================\n";
    
    return 0;
}