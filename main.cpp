#include <iostream>
#include <vector>

using vvi = std::vector<std::vector<double>>;
using vi = std::vector<double>;

bool check_max_dg(int n, vvi &matrix)  {// TODO: wtf как поравнять и в строках и в столбцах
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i != j && matrix[i][i] <= matrix[i][j]) {
                return false;
            }
        }
    }
    return true;
}

int splitLU(int n, vvi &matrix, vvi &L, vvi &U) {
    for (int j = 0; j < n; ++j) {
        U[1][j] = matrix[1][j];
    }
    for (int j = 1; j < n; ++j) {
        L[j][1] = matrix[j][1] / U[1][j];
    }
    for (int i = 0; i < n; ++i){
        for (int j = i; j < n; ++j) {

        }
    }
}

int main() {
    int n;
    double elem;
    std::cout << "Input number of variables:\n";
    std::cin >> n;
    std::cout << "Input matrix:\n";
    vvi matrixA(n, (vi(n)));
    vi b(n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            std::cin >> elem;
            matrixA[i][j] = elem;
        }
        std::cin >> elem;
        b[i] = elem;
    }
//    check_max_dg(n, matrixA);
    //TODO: функция по разложению на матрицы LU
}