#ifndef LINEARALGEBRA_LIBRARY_H
#define LINEARALGEBRA_LIBRARY_H

#include <algorithm>
#include <cstdio>
#include <iostream>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace la {
class matrix {
    std::vector<double> data_;
    size_t amount_of_rows_;
    size_t amount_of_columns_;

public:
    matrix(size_t const &rows, size_t const &columns)
        : data_(rows * columns), amount_of_rows_(rows), amount_of_columns_(columns) {}

    matrix(matrix const &rhs) = default;

    matrix() = default;

    double const &operator()(size_t const &row, size_t const &column) const;

    double &operator()(size_t const &row, size_t const &column);

    matrix &operator+=(matrix const &rhs);

    friend matrix operator+(matrix const &lhs, matrix const &rhs);

    matrix &operator-=(matrix const &rhs);

    friend matrix operator-(matrix const &lhs, matrix const &rhs);

    matrix &operator/=(double const &rhs);

    friend matrix operator/(matrix const &lhs, double const &rhs);

    matrix &operator*=(double const &rhs);

    friend matrix operator*(matrix const &lhs, double const &rhs);

    friend matrix operator*(double const &lhs, matrix const &rhs);

    matrix &operator*=(matrix const &rhs);

    friend matrix operator*(matrix const &lhs, matrix &rhs);

    friend bool operator==(matrix const &lhs, matrix const &rhs);

    friend bool operator!=(matrix const &lhs, matrix const &rhs);

    double determinant() const;

    matrix &transpose();

    matrix &inverse();

    size_t height() const;

    size_t length() const;

    friend std::ostream &operator<<(std::ostream &cout, matrix const &a);

    friend std::istream &operator>>(std::istream &cin, matrix &a);
};

matrix transpose(matrix const &a);

matrix inverse(matrix const &a);

matrix E(size_t i);

class tensor {
    std::vector<double> data_;
    size_t amount_of_p_;
    size_t amount_of_q_;
    size_t dimension_;

    size_t binary_pow(size_t const &dim, size_t const &power);

    double sum_el(std::vector<size_t> const &old, std::vector<size_t> const &real, matrix T_matrix, matrix S_matrix);

    bool change_round(std::vector<size_t> &round) const;

public:
    tensor(size_t const &dimension, size_t const &p, size_t const &q)
        : data_(binary_pow(dimension, (p + q))), amount_of_p_(p), amount_of_q_(q), dimension_(dimension) {}

    tensor(tensor const &rhs) = default;

    tensor() = default;

    double const &operator()(std::vector<size_t> const &v) const;

    double &operator()(std::vector<size_t> const &v);

    double const &operator()() const;

    double &operator()();

    tensor &contraction(size_t p_num, size_t q_num);

    tensor &change_basis(matrix const &T_matrix);

    tensor &sym(std::vector<bool> const &v);

    tensor &asym(std::vector<bool> const &v);

    friend std::ostream &operator<<(std::ostream &cout, tensor const &a);

    friend std::istream &operator>>(std::istream &cin, tensor &a);
};

tensor contraction(tensor const &a, size_t p_num, size_t q_num);

tensor change_basis(tensor const &a, matrix const &T_matrix);

tensor sym(tensor const &a, std::vector<bool> const &v);

tensor asym(tensor const &a, std::vector<bool> const &v);
}  // namespace la
#endif