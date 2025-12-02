#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <functional>

using namespace std;

struct Solution {
    vector<double> x;
    vector<double> y;
    vector<double> dy;
};

struct EulerCauchySolution {
    vector<double> x;
    vector<double> y;
    vector<double> y_pred;  // прогнозные значения
    vector<double> delta_y; // приращения
};

struct ImprovedEulerSolution {
    vector<double> x;
    vector<double> y;
    vector<double> y_half;  // значения в средней точке
    vector<double> delta_y; // приращения
};

// Точное решение: y = (e^(x²) + e^(-x²) - 1)e^(x²)
double exact_solution(double x) {
    double x_squared = x * x;
    double exp_x2 = exp(x_squared);
    double exp_minus_x2 = exp(-x_squared);
    return (exp_x2 + exp_minus_x2 - 1) * exp_x2;
}

// Правая часть уравнения: y" = A*x*y' - (4x² - 3)*y + e^(x²)
// Уравнение: y" - A*x*y' + (4x² - 3)y = e^(x²)
// Где A = 2 (вероятно, из контекста)
double f(double x, double y, double z) {
    // Здесь z = y'
    double A = 2.0; // Из уравнения видно A=2
    double x_squared = x * x;
    return A * x * z - (4 * x_squared - 3) * y + exp(x_squared);
}

// Для методов Эйлера, Эйлера-Коши и улучшенного Эйлера нужно f(x, y)
// Преобразуем уравнение второго порядка к системе первого порядка
// y' = z
// z' = f(x, y, z)

// ================== МЕТОД ЭЙЛЕРА (явный) ДЛЯ СИСТЕМЫ ==================
Solution euler_method_system(double x0, double y0, double z0, double h, int n) {
    Solution sol;
    sol.x.resize(n + 1);
    sol.y.resize(n + 1);
    sol.dy.resize(n + 1); // здесь храним z = y'

    sol.x[0] = x0;
    sol.y[0] = y0;
    sol.dy[0] = z0;

    for (int k = 0; k < n; k++) {
        double xk = sol.x[k];
        double yk = sol.y[k];
        double zk = sol.dy[k];

        sol.x[k + 1] = xk + h;
        sol.y[k + 1] = yk + h * zk; // y' = z
        sol.dy[k + 1] = zk + h * f(xk, yk, zk); // z' = f(x, y, z)
    }

    return sol;
}

// ================== МЕТОД ЭЙЛЕРА-КОШИ ДЛЯ СИСТЕМЫ ==================
EulerCauchySolution euler_cauchy_method_system(double x0, double y0, double z0, double h, int n) {
    EulerCauchySolution sol;
    sol.x.resize(n + 1);
    sol.y.resize(n + 1);
    sol.y_pred.resize(n + 1);
    sol.delta_y.resize(n + 1);

    // Для системы храним также производные
    vector<double> z(n + 1);
    vector<double> z_pred(n + 1);

    sol.x[0] = x0;
    sol.y[0] = y0;
    z[0] = z0;

    for (int k = 0; k < n; k++) {
        double xk = sol.x[k];
        double yk = sol.y[k];
        double zk = z[k];

        // Прогноз (этап 1) - метод Эйлера
        double y_pred = yk + h * zk;
        double z_pred_val = zk + h * f(xk, yk, zk);

        sol.y_pred[k + 1] = y_pred;
        z_pred[k + 1] = z_pred_val;

        // Коррекция (этап 2) - усреднение
        double x_next = xk + h;
        double y_next = yk + 0.5 * h * (zk + z_pred_val);
        double z_next = zk + 0.5 * h * (f(xk, yk, zk) + f(x_next, y_pred, z_pred_val));

        sol.x[k + 1] = x_next;
        sol.y[k + 1] = y_next;
        sol.delta_y[k] = y_next - yk;
        z[k + 1] = z_next;
    }

    return sol;
}

// ================== ПЕРВЫЙ УЛУЧШЕННЫЙ МЕТОД ЭЙЛЕРА ДЛЯ СИСТЕМЫ ==================
ImprovedEulerSolution improved_euler_method_system(double x0, double y0, double z0, double h, int n) {
    ImprovedEulerSolution sol;
    sol.x.resize(n + 1);
    sol.y.resize(n + 1);
    sol.y_half.resize(n + 1);
    sol.delta_y.resize(n + 1);

    // Для системы храним также производные
    vector<double> z(n + 1);
    vector<double> z_half(n + 1);

    sol.x[0] = x0;
    sol.y[0] = y0;
    z[0] = z0;

    for (int k = 0; k < n; k++) {
        double xk = sol.x[k];
        double yk = sol.y[k];
        double zk = z[k];

        // Шаг на половинном интервале
        double x_half = xk + h/2;
        double y_half = yk + (h/2) * zk;
        double z_half_val = zk + (h/2) * f(xk, yk, zk);

        // Полный шаг с использованием производных в средней точке
        double y_next = yk + h * z_half_val;
        double z_next = zk + h * f(x_half, y_half, z_half_val);

        sol.x[k + 1] = xk + h;
        sol.y[k + 1] = y_next;
        sol.y_half[k] = y_half;
        sol.delta_y[k] = y_next - yk;
        z[k + 1] = z_next;
    }

    return sol;
}

// ================== МЕТОД РУНГЕ-КУТТЫ 4-го порядка ДЛЯ СИСТЕМЫ ==================
Solution runge_kutta_4_method_system(double x0, double y0, double z0, double h, int n) {
    Solution sol;
    sol.x.resize(n + 1);
    sol.y.resize(n + 1);
    sol.dy.resize(n + 1); // здесь храним z = y'

    sol.x[0] = x0;
    sol.y[0] = y0;
    sol.dy[0] = z0;

    for (int k = 0; k < n; k++) {
        double xk = sol.x[k];
        double yk = sol.y[k];
        double zk = sol.dy[k];

        double K1_y = h * zk;
        double K1_z = h * f(xk, yk, zk);

        double K2_y = h * (zk + K1_z/2);
        double K2_z = h * f(xk + h/2, yk + K1_y/2, zk + K1_z/2);

        double K3_y = h * (zk + K2_z/2);
        double K3_z = h * f(xk + h/2, yk + K2_y/2, zk + K2_z/2);

        double K4_y = h * (zk + K3_z);
        double K4_z = h * f(xk + h, yk + K3_y, zk + K3_z);

        sol.x[k + 1] = xk + h;
        sol.y[k + 1] = yk + (K1_y + 2*K2_y + 2*K3_y + K4_y) / 6;
        sol.dy[k + 1] = zk + (K1_z + 2*K2_z + 2*K3_z + K4_z) / 6;
    }

    return sol;
}

// ================== МЕТОД АДАМСА 4-го порядка ДЛЯ СИСТЕМЫ ==================
Solution adams_4_method_system(double x0, double y0, double z0, double h, int n) {
    Solution sol;
    sol.x.resize(n + 1);
    sol.y.resize(n + 1);
    sol.dy.resize(n + 1); // здесь храним z = y' и f(x,y,z)
    vector<double> f_vals(n + 1); // храним f(x,y,z)

    // Запускаем метод Рунге-Кутты для получения первых 4 точек
    Solution rk_start = runge_kutta_4_method_system(x0, y0, z0, h, 3);

    for (int i = 0; i < 4; i++) {
        sol.x[i] = rk_start.x[i];
        sol.y[i] = rk_start.y[i];
        sol.dy[i] = rk_start.dy[i];
        f_vals[i] = f(sol.x[i], sol.y[i], sol.dy[i]);
    }

    // Продолжаем методом Адамса, начиная с 4-й точки
    for (int k = 3; k < n; k++) {
        sol.x[k + 1] = sol.x[k] + h;

        // Формула Адамса для y
        sol.y[k + 1] = sol.y[k] + h/24.0 * (55*sol.dy[k] - 59*sol.dy[k-1]
                                          + 37*sol.dy[k-2] - 9*sol.dy[k-3]);

        // Сначала используем прогноз для z
        double z_pred = sol.dy[k] + h/24.0 * (55*f_vals[k] - 59*f_vals[k-1]
                                           + 37*f_vals[k-2] - 9*f_vals[k-3]);

        // Затем вычисляем f в новой точке
        f_vals[k + 1] = f(sol.x[k + 1], sol.y[k + 1], z_pred);

        // Коррекция для z
        sol.dy[k + 1] = sol.dy[k] + h/24.0 * (55*f_vals[k] - 59*f_vals[k-1]
                                           + 37*f_vals[k-2] - 9*f_vals[k-3]);

        // Уточняем f с новым значением z
        f_vals[k + 1] = f(sol.x[k + 1], sol.y[k + 1], sol.dy[k + 1]);
    }

    return sol;
}

// ================== ВЫВОД ТАБЛИЦ ==================
void print_simple_table(const Solution& sol, const string& method_name) {
    cout << "\nTable by method of " << method_name << "\n";
    cout << "==============================================================\n";
    cout << setw(5) << "k" << setw(12) << "x_k" << setw(15) << "y_k"
         << setw(15) << "y_actual" << setw(15) << "error" << endl;
    cout << "==============================================================\n";

    for (size_t k = 0; k < sol.x.size(); k++) {
        double y_exact = exact_solution(sol.x[k]);
        double error = abs(y_exact - sol.y[k]);

        cout << setw(5) << k
             << setw(12) << fixed << setprecision(6) << sol.x[k]
             << setw(15) << sol.y[k]
             << setw(15) << y_exact
             << setw(15) << setprecision(6) << error << endl;
    }
    cout << "==============================================================\n";
}

void print_euler_cauchy_table(const EulerCauchySolution& sol, const string& method_name) {
    cout << "\nTable by method of " << method_name << "\n";
    cout << "============================================================================================\n";
    cout << setw(5) << "k" << setw(12) << "x_k" << setw(15) << "y_k"
         << setw(15) << "y_pred" << setw(15) << "Δy_k"
         << setw(15) << "y_actual" << setw(15) << "error" << endl;
    cout << "============================================================================================\n";

    for (size_t k = 0; k < sol.x.size(); k++) {
        double y_exact = exact_solution(sol.x[k]);
        double error = abs(y_exact - sol.y[k]);

        cout << setw(5) << k
             << setw(12) << fixed << setprecision(6) << sol.x[k]
             << setw(15) << sol.y[k];

        if (k < sol.y_pred.size() - 1) {
            cout << setw(15) << sol.y_pred[k + 1]
                 << setw(15) << sol.delta_y[k];
        } else {
            cout << setw(15) << "-"
                 << setw(15) << "-";
        }

        cout << setw(15) << y_exact
             << setw(15) << setprecision(6) << error << endl;
    }
    cout << "============================================================================================\n";
}

void print_improved_euler_table(const ImprovedEulerSolution& sol, const string& method_name) {
    cout << "\nTable by method of " << method_name << "\n";
    cout << "============================================================================================\n";
    cout << setw(5) << "k" << setw(12) << "x_k" << setw(15) << "y_k"
         << setw(15) << "y_{k+1/2}" << setw(15) << "Δy_k"
         << setw(15) << "y_actual" << setw(15) << "error" << endl;
    cout << "============================================================================================\n";

    for (size_t k = 0; k < sol.x.size(); k++) {
        double y_exact = exact_solution(sol.x[k]);
        double error = abs(y_exact - sol.y[k]);

        cout << setw(5) << k
             << setw(12) << fixed << setprecision(6) << sol.x[k]
             << setw(15) << sol.y[k];

        if (k < sol.y_half.size()) {
            cout << setw(15) << sol.y_half[k]
                 << setw(15) << sol.delta_y[k];
        } else {
            cout << setw(15) << "-"
                 << setw(15) << "-";
        }

        cout << setw(15) << y_exact
             << setw(15) << setprecision(6) << error << endl;
    }
    cout << "============================================================================================\n";
}

// ================== ОЦЕНКА ПОГРЕШНОСТИ ==================
double runge_romberg_error(double y_h, double y_2h, int p) {
    return abs(y_h - y_2h) / (pow(2, p) - 1);
}

// ================== ОСНОВНАЯ ПРОГРАММА ==================
int main() {
    // Параметры задачи из фото
    double x0 = 0.0;
    double y0 = 1.0;       // y(0) = 1
    double z0 = 0.0;       // y'(0) = 0
    double b = 1.0;        // x ∈ [0, 1]
    double h = 0.1;        // шаг
    int n = (int)((b - x0) / h);

    cout << "==============================================\n";
    cout << "   SOLUTION OF CAUCHY PROBLEM\n";
    cout << "   y'' - 2xy' + (4x² - 3)y = e^(x²)\n";
    cout << "   y(0) = 1, y'(0) = 0\n";
    cout << "   on interval [0, 1.0]\n";
    cout << "   step h = " << h << "\n";
    cout << "   Exact solution: y(x) = (e^(x²) + e^(-x²) - 1)e^(x²)\n";
    cout << "==============================================\n";

    // 1. Метод Эйлера (явный) для системы
    cout << "\n1. EXPLICIT EULER METHOD FOR SYSTEM\n";
    Solution sol_euler = euler_method_system(x0, y0, z0, h, n);
    print_simple_table(sol_euler, "explicit Euler");

    // 2. Метод Эйлера-Коши для системы
    cout << "\n2. EULER-CAUCHY METHOD FOR SYSTEM\n";
    EulerCauchySolution sol_euler_cauchy = euler_cauchy_method_system(x0, y0, z0, h, n);
    print_euler_cauchy_table(sol_euler_cauchy, "Euler-Cauchy");

    // 3. Первый улучшенный метод Эйлера для системы
    cout << "\n3. FIRST IMPROVED EULER METHOD FOR SYSTEM\n";
    ImprovedEulerSolution sol_improved_euler = improved_euler_method_system(x0, y0, z0, h, n);
    print_improved_euler_table(sol_improved_euler, "first improved Euler");

    // 4. Метод Рунге-Кутты 4-го порядка для системы
    cout << "\n4. RUNGE-KUTTA 4th ORDER METHOD FOR SYSTEM\n";
    Solution sol_rk4 = runge_kutta_4_method_system(x0, y0, z0, h, n);
    print_simple_table(sol_rk4, "Runge-Kutta 4th order");

    // 5. Метод Адамса 4-го порядка для системы
    cout << "\n5. ADAMS 4th ORDER METHOD FOR SYSTEM\n";
    Solution sol_adams = adams_4_method_system(x0, y0, z0, h, n);
    print_simple_table(sol_adams, "Adams 4th order");

    // ================== ОЦЕНКА ПОГРЕШНОСТИ ==================
    cout << "\n6. ERROR ESTIMATION\n";
    cout << "==============================================\n";

    double h2 = h * 2;
    int n2 = (int)((b - x0) / h2);

    // Для метода Эйлера
    Solution sol_euler_h = euler_method_system(x0, y0, z0, h, n);
    Solution sol_euler_2h = euler_method_system(x0, y0, z0, h2, n2);
    double error_euler = runge_romberg_error(sol_euler_h.y[n], sol_euler_2h.y[n2], 1);

    // Для метода Эйлера-Коши
    EulerCauchySolution sol_ec_h = euler_cauchy_method_system(x0, y0, z0, h, n);
    EulerCauchySolution sol_ec_2h = euler_cauchy_method_system(x0, y0, z0, h2, n2);
    double error_ec = runge_romberg_error(sol_ec_h.y[n], sol_ec_2h.y[n2], 2);

    // Для улучшенного метода Эйлера
    ImprovedEulerSolution sol_ie_h = improved_euler_method_system(x0, y0, z0, h, n);
    ImprovedEulerSolution sol_ie_2h = improved_euler_method_system(x0, y0, z0, h2, n2);
    double error_ie = runge_romberg_error(sol_ie_h.y[n], sol_ie_2h.y[n2], 2);

    // Для метода Рунге-Кутты
    Solution sol_rk4_h = runge_kutta_4_method_system(x0, y0, z0, h, n);
    Solution sol_rk4_2h = runge_kutta_4_method_system(x0, y0, z0, h2, n2);
    double error_rk4 = runge_romberg_error(sol_rk4_h.y[n], sol_rk4_2h.y[n2], 4);

    cout << "Error estimation by Runge-Romberg method:\n";
    cout << "------------------------------------------\n";
    cout << "Euler method (p=1):           " << setprecision(6) << error_euler << endl;
    cout << "Euler-Cauchy method (p=2):    " << error_ec << endl;
    cout << "Improved Euler method (p=2):  " << error_ie << endl;
    cout << "Runge-Kutta 4 (p=4):          " << error_rk4 << endl;

    // ================== СРАВНЕНИЕ С ТОЧНЫМ РЕШЕНИЕМ ==================
    cout << "\n7. COMPARISON WITH EXACT SOLUTION\n";
    cout << "==============================================\n";

    double y_exact_end = exact_solution(b);

    cout << "Exact solution y(" << b << ") = " << fixed << setprecision(9) << y_exact_end << endl;
    cout << "\nComparison table at x = " << b << ":\n";
    cout << "=====================================================================================\n";
    cout << setw(25) << "Method" << setw(8) << "y(" << b << ")"
         << setw(20) << "Absolute error" << setw(20) << "RR error estimate" << endl;
    cout << "=====================================================================================\n";

    vector<pair<string, double>> results = {
        {"Explicit Euler", sol_euler.y[n]},
        {"Euler-Cauchy", sol_euler_cauchy.y[n]},
        {"Improved Euler", sol_improved_euler.y[n]},
        {"Runge-Kutta 4", sol_rk4.y[n]},
        {"Adams 4", sol_adams.y[n]}
    };

    vector<double> rr_errors = {error_euler, error_ec, error_ie, error_rk4, 0};

    for (size_t i = 0; i < results.size(); i++) {
        double y_val = results[i].second;
        double abs_error = abs(y_exact_end - y_val);
        double rr_error = (i < rr_errors.size()) ? rr_errors[i] : 0;

        cout << setw(25) << results[i].first
             << setw(20) << fixed << setprecision(9) << y_val
             << setw(20) << setprecision(6) << abs_error
             << setw(20) << rr_error << endl;
    }
    cout << "=====================================================================================\n";

    // ================== КОМПАКТНАЯ СВОДКА (как в методичке) ==================
    cout << "\n8. COMPACT SUMMARY TABLE\n";
    cout << "==========================================================================\n";
    cout << "k" << setw(12) << "x_k" << setw(12) << "Euler" << setw(12) << "Eul-Cauchy"
         << setw(12) << "Imp.Euler" << setw(12) << "RK4" << setw(12) << "Exact" << endl;
    cout << "==========================================================================\n";

    for (int k = 0; k <= n; k++) {
        double xk = x0 + k * h;
        cout << fixed << setprecision(1) << setw(1) << k
             << setprecision(5) << setw(12) << xk
             << setw(12) << sol_euler.y[k]
             << setw(12) << sol_euler_cauchy.y[k]
             << setw(12) << sol_improved_euler.y[k]
             << setw(12) << sol_rk4.y[k]
             << setw(12) << exact_solution(xk) << endl;
    }
    cout << "==========================================================================\n";

    // ================== ТАБЛИЦА ПРОИЗВОДНЫХ ==================
    cout << "\n9. DERIVATIVES TABLE (y' values)\n";
    cout << "=====================================================\n";
    cout << setw(8) << "x" << setw(15) << "Euler y'" << setw(15) << "RK4 y'"
         << setw(15) << "Adams y'" << endl;
    cout << "=====================================================\n";

    for (int k = 0; k <= n; k += 2) { // Выводим через точку
        cout << fixed << setprecision(2) << setw(8) << sol_euler.x[k]
             << setprecision(8) << setw(15) << sol_euler.dy[k]
             << setw(15) << sol_rk4.dy[k]
             << setw(15) << sol_adams.dy[k] << endl;
    }
    cout << "=====================================================\n";
    
    return 0;
}