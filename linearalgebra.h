#ifndef LINEARALGEBRA_LIBRARY_H
#define LINEARALGEBRA_LIBRARY_H

#include <vector>
#include <algorithm>
#include <cstdio>
#include <stdexcept>
#include <type_traits>

namespace la {
    class matrix {
        std::vector<double> data_;
        size_t amount_of_rows_;
        size_t amount_of_columns_;
     public:
        matrix(size_t const &rows, size_t const &columns) : data_(rows * columns),
                                                            amount_of_rows_(rows),
                                                            amount_of_columns_(columns) {}

        matrix(matrix const &rhs) = default;

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

        double sum_el(std::vector<size_t> const &old, std::vector<size_t> const &real,
                      matrix T_matrix, matrix S_matrix);

        bool change_round(std::vector<size_t> &round);

     public:
        tensor(size_t const &dimension, size_t const &p, size_t const &q) : data_(binary_pow(dimension, (p + q))),
                                                                            amount_of_p_(p),
                                                                            amount_of_q_(q),
                                                                            dimension_(dimension) {}

        tensor(tensor const &rhs) = default;

        double const &operator()(std::vector<size_t> const &v) const;

        double &operator()(std::vector<size_t> const &v);

        tensor &contraction(size_t p_num, size_t q_num);

        tensor change_basis(matrix const &T_matrix);
    };

    tensor contraction(tensor const &a, size_t p_num, size_t q_num);
}
#endif