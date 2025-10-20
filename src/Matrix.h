#pragma once
#include <bits/stdc++.h>

template<class T>
class matrix {
    std::vector<std::vector<T>> data;
    size_t rows;
    size_t cols;

    class RowProxy {
        std::vector<T>& row;
    public:
        explicit RowProxy(std::vector<T>& r);
        T& operator[](size_t col);
    };

    class ConstRowProxy {
        const std::vector<T>& row;
    public:
        explicit ConstRowProxy(const std::vector<T>& r);
        const T& operator[](size_t col) const;
    };

public:
    matrix();
    matrix(size_t n, size_t m);
    matrix(size_t n, size_t m, T val);
    matrix(std::initializer_list<std::initializer_list<T>> init);
    matrix(const matrix& other);
    ~matrix();

    matrix& operator=(const matrix& other);
    matrix& operator=(std::initializer_list<std::initializer_list<T>> init);
    matrix& operator+=(const matrix& other);
    matrix& operator-=(const matrix& other);
    matrix& operator*=(const T& other);
    matrix& operator/=(const T& other);
    matrix& operator*=(const matrix& other);

    matrix operator*(const matrix& other) const;
    std::vector<T> operator*(const std::vector<T>& other) const;
    matrix operator+(const matrix& other) const;
    matrix operator-(const matrix& other) const;
    matrix operator*(const T& other) const;
    matrix operator/(const T& other) const;

    T& operator()(size_t row, size_t col);
    RowProxy operator[](size_t row);
    ConstRowProxy operator[](size_t row) const;
    bool operator==(const matrix& other);

    static matrix<T> setE(const int& n);
    matrix transpose() const;
    matrix inverse() const;
    matrix& resize(size_t new_rows, size_t new_cols);

    size_t getRows() const;
    size_t getCols() const;
};


template<class T>
matrix<T>::matrix() : data(), rows(0), cols(0) {}

template<class T>
matrix<T>::matrix(const size_t n, const size_t m) : data(n, std::vector<T> (m, 0)), rows(n), cols(m) {}

template<class T>
matrix<T>::matrix(size_t n, size_t m, T val) : data(n, vector<T> (m, val)), rows(n), cols(m) {}

template<class T>
matrix<T>::matrix(std::initializer_list<std::initializer_list<T>> init) : rows(init.size()), cols(0) {

    for (const auto& row_list : init) {
        if (row_list.size() > cols) {
            cols = row_list.size();
        }
    }

    data.resize(rows);
    size_t i = 0;
    for (const auto& row_list : init) {
        data[i].reserve(cols);

        for (const auto& elem : row_list) {
            data[i].push_back(elem);
        }
        while (data[i].size() < cols) {
            data[i].push_back(T{});
        }
        ++i;
    }
}

template<class T>
matrix<T>::matrix(const matrix &other) : data(other.data), rows(other.rows), cols(other.cols) {}

template<class T>
matrix<T>::~matrix() = default;

template<class T>
matrix<T>::RowProxy::RowProxy (std::vector<T>& r) : row(r) {}
template<class T>
T& matrix<T>::RowProxy::operator[] (const size_t col) {
    if (col >= row.size()) {
        throw std::out_of_range("matrix::RowProxy::operator[]");
    }
    return row[col];
}

template<class T>
matrix<T>::ConstRowProxy::ConstRowProxy (const std::vector<T>& r) : row(r) {}
template<class T>
const T& matrix<T>::ConstRowProxy::operator[](size_t col) const {
    if (col >= row.size()) {
        throw std::out_of_range("matrix::ConstRowProxy::operator[]");
    }
    return row[col];
}

template<class T>
matrix<T>& matrix<T>::operator=(const matrix &other) {
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            data[i][j] = other.data[i][j];
        }
        for (size_t j = cols; j < other.cols; j++) {
            data[i].push_back(other.data[i][j]);
        }
        if (cols > other.cols) {
            data[i].erase(data[i].begin() + other.cols, data[i].end());
        }
    }

    for (size_t i = rows; i < other.rows; i++) {
        data.push_back(other.data[i]);
    }
    if (rows > other.rows) {
        data.erase(data.begin() + other.rows, data.end());
    }

    rows = other.rows;
    cols = other.cols;
    return *this;
}

template<class T>
matrix<T>& matrix<T>::operator=(std::initializer_list<std::initializer_list<T>> init) {
    size_t new_rows = init.size();
    size_t new_cols = 0;

    for (const auto& row_list : init) {
        if (row_list.size() > new_cols) {
            new_cols = row_list.size();
        }
    }

    data.resize(new_rows);
    size_t i = 0;
    for (const auto& row_list : init) {
        data[i].clear();
        for (const auto& elem : row_list) {
            data[i].push_back(elem);
        }
        while (data[i].size() < new_cols) {
            data[i].push_back(T{});
        }
        ++i;
    }

    rows = new_rows;
    cols = new_cols;
    return *this;
}

template<class T>
matrix<T>& matrix<T>::operator+=(const matrix &other) {
    if (cols != other.cols || rows != other.rows) {
        throw std::invalid_argument("matrix size mismatch in +=");
    }
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            data[i][j] += other.data[i][j];
        }
    }
    return *this;
}

template<class T>
matrix<T>& matrix<T>::operator-=(const matrix &other) {
    if (cols != other.cols || rows != other.rows) {
        throw std::invalid_argument("matrix size mismatch in -=");
    }
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            data[i][j] -= other.data[i][j];
        }
    }
    return *this;
}

template<class T>
matrix<T>& matrix<T>::operator*=(const T &other) {
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            data[i][j] *= other;
        }
    }
    return *this;
}

template<class T>
matrix<T>& matrix<T>::operator/=(const T &other) {
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            data[i][j] /= other;
        }
    }
    return *this;
}

template<class T>
matrix<T> matrix<T>::operator*(const matrix &other) const {
    if (cols != other.rows) {
        throw std::invalid_argument("matrix size mismatch in *=");
    }
    matrix tmp(rows, other.cols);

    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < other.cols; j++) {
            for (size_t k = 0; k < cols; k++) {
                tmp[i][j] += data[i][k] * other.data[k][j];
            }
        }
    }

    return tmp;
}

template<class T>
std::vector<T> matrix<T>::operator*(const std::vector<T> &other) const {
    if (cols != other.size()) {
        throw std::invalid_argument("matrix size mismatch in *=");
    }
    std::vector<T> tmp(rows);

    for (size_t i = 0; i < rows; i++) {
        for (size_t k = 0; k < cols; k++) {
            tmp[i] += data[i][k] * other[k];
        }
    }

    return tmp;
}

template<class T>
matrix<T>& matrix<T>::operator*=(const matrix &other) {
    *this = *this * other;
    return *this;
}

template<class T>
matrix<T> matrix<T>::operator+(const matrix &other) const {
    matrix tmp(*this);
    tmp += other;
    return tmp;
}

template<class T>
matrix<T> matrix<T>::operator-(const matrix &other) const {
    matrix tmp(*this);
    tmp -= other;
    return tmp;
}

template<class T>
matrix<T> matrix<T>::operator*(const T &other) const {
    matrix tmp(*this);
    tmp *= other;
    return tmp;
}

template<class T>
matrix<T> matrix<T>::operator/(const T &other) const {
    matrix tmp(*this);
    tmp *= other;
    return tmp;
}

template<class T>
T& matrix<T>::operator()(size_t row, size_t col) {
    if (row >= rows || col >= cols) {
        throw std::invalid_argument("matrix size mismatch in ()");
    }
    return data[row][col];
}

template<class T>
bool matrix<T>::operator==(const matrix &other) {
    if (cols != other.cols || rows != other.rows) {
        return false;
    }
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            if (data[i][j] != other.data[i][j]) {
                return false;
            }
        }
    }
    return true;
}

template<class T>
typename matrix<T>::RowProxy matrix<T>::operator[](size_t row) {
    if (row >= rows) {
        throw std::invalid_argument("matrix size mismatch in []");
    }
    return RowProxy(data[row]);
}
template<class T>
typename matrix<T>::ConstRowProxy matrix<T>::operator[](size_t row) const {
    if (row >= rows) {
        throw std::out_of_range("matrix row index out of range");
    }
    return ConstRowProxy(data[row]);
}

template<class T>
matrix<T> matrix<T>::setE(const int& n) {
    matrix<T> tmp(n);
    for (int i = 0; i < n; i++) {
        tmp[i][i] = 1;
    }
    return tmp;
}

template<class T>
matrix<T> matrix<T>::transpose() const {
//    if (cols != rows) {
//        throw std::invalid_argument("matrix size mismatch in transpose()");
//    }
    matrix tmp(cols, rows);
    for (size_t i = 0; i < cols; i++) {
        for (size_t j = 0; j < rows; j++) {
            tmp.data[i][j] = data[j][i];
        }
    }
    return tmp;
}

template<class T>
matrix<T> matrix<T>::inverse() const {
    if (cols != rows) {
        throw std::invalid_argument("matrix size mismatch in ()");
    }
    if (cols != 2) {
        throw std::invalid_argument("so difficult for invers method");
    }
    matrix tmp(cols, cols);

    double det = data[0][0] * data[1][1] - data[0][1] * data[1][0];
    tmp.data[0][0] = data[1][1] / det;
    tmp.data[0][1] = -data[0][1] / det;
    tmp.data[1][0] = -data[1][0] / det;
    tmp.data[1][1] = data[0][0] / det;

    return tmp;
}

template<class T>
matrix<T> &matrix<T>::resize(size_t new_rows, size_t new_cols) {
    if (new_rows > rows) {
        data.insert(data.end(), new_rows - rows, std::vector<T>(new_cols, 0));
    }
    if (new_rows < rows) {
        data.erase(data.begin() + new_rows, data.end());
    }

    for (size_t i = 0; i < rows; i++) {
        if (new_cols > cols) {
            data[i].insert(data[i].end(), new_cols - cols, 0);
        }
        if (new_cols < cols) {
            data[i].erase(data[i].begin() + new_cols, data[i].end());
        }
    }
    return *this;
}

template<class T>
size_t matrix<T>::getRows() const { return rows; }

template<class T>
size_t matrix<T>::getCols() const { return cols; }

template<class T>
std::ostream& operator<<(std::ostream& os, const matrix<T>& matrix) {
    for (size_t i = 0; i < matrix.getRows(); i++) {
        for (size_t j = 0; j < matrix.getCols(); j++) {
            if (fabs(matrix[i][j]) < 0.00001) {
                os << std::fixed << std::setprecision(2) << 0.0;
            } else {
                os << std::fixed << std::setprecision(2) << matrix[i][j];
            }

            if (j < matrix.getCols() - 1) {
                os << " ";
            }
        }
        if (i < matrix.getRows() - 1) {
            os << "\n";
        }
    }
    return os;
}