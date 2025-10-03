#include <iostream>
#include <vector>

template <class T>
class matrix {
    size_t rows;
    size_t cols;
    std::vector<std::vector<T>> data;

    class RowProxy {
    private:
        std::vector<int>& row;
    public:
        explicit RowProxy(std::vector<int>& r);
        int& operator[](size_t col);
    };

public:
    matrix();

    matrix(size_t n, size_t m);

    matrix(size_t n, size_t m, T val);

    matrix(std::initializer_list<T> init);

    matrix(const matrix &other);

    ~matrix();

    matrix& operator=(const matrix &other);

    matrix& operator=(std::initializer_list<T> init);

    matrix& operator+=(matrix &other);

    matrix& operator-=(matrix &other);

    matrix& operator*=(matrix &other);

    std::vector<T>& operator*=(std::vector<T> &other);

    matrix operator+(matrix &other);

    matrix operator-(matrix &other);

    matrix operator*(matrix &other);

    std::vector<T> operator*(std::vector<T> &other);

    bool operator==(matrix &other);

    T& operator()(size_t row, size_t col);

    RowProxy& operator[](size_t row);
};