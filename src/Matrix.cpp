#include "Matrix.h"

template<class T>
matrix<T>::matrix() : data(), rows(0), cols(0) {}

template<class T>
matrix<T>::matrix(const size_t n, const size_t m) : data(n, vector<T> (m, 0)), rows(n), cols(m) {}

template<class T>
matrix<T>::matrix(size_t n, size_t m, T val) : data(n, vector<T> (m, val)), rows(n), cols(m) {}

template<class T>
matrix<T>::matrix(std::initializer_list<T> init) : data(init), rows(init.size()), cols(0) {
    for(int i = 0; i < rows; i++) {
        if (cols < data[i].size()) {
            cols = data[i].size();
        }
    }
    for (int i = 0; i < rows; i++) {
        size_t diff;
        diff = cols - data[i].size();
        data[i].insert(data[i].end(), diff, 0);
    }
}

template<class T>
matrix<T>::matrix(const matrix &other) : data(other.data), rows(other.rows), cols(other.cols) {}

template<class T>
matrix<T>::~matrix() {
    delete [] data;
}

template<class T>
matrix<T>& matrix<T>::operator=(const matrix<T> &other) {
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
matrix<T>& matrix<T>::operator=(std::initializer_list<T> init) {
    size_t other_rows = init.size();
    size_t other_cols = init[0].size();
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            data[i][j] = init[i][j];
        }
        for (size_t j = cols; j < other_cols; j++) {
            data[i].push_back(init[i][j]);
        }
        if (cols > other_cols) {
            data[i].erase(data[i].begin() + other_cols, data[i].end());
        }
    }
    for (size_t i = rows; i < other_rows; i++) {
        data.push_back(init[i]);
    }
    if (rows > other_rows) {
        data.erase(data.begin() + other_rows, data.end());
    }
    rows = other_rows;
    cols = other_cols;
    return *this;
}



