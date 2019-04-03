#include "linearalgebra.h"

namespace la {
    double const &matrix::operator()(size_t const &row, size_t const &column) const {
        if (row >= amount_of_rows_ || column >= amount_of_columns_) {
            char error[256];
            std::sprintf(error,
                         "Adress to [%zu][%zu], size of matrix is [%zu][%zu]",
                         row, column, amount_of_rows_, amount_of_columns_);
            throw std::out_of_range(error);
        }
        return data_[row * amount_of_columns_ + column];
    }

    double &matrix::operator()(size_t const &row, size_t const &column) {
        return const_cast<double &>(const_cast<const matrix *>(this)->operator()(row, column));
    }

    matrix &matrix::operator+=(matrix const &rhs) {
        if (amount_of_rows_ != rhs.amount_of_rows_ || amount_of_columns_ != rhs.amount_of_columns_) {
            char error[256];
            std::sprintf(error,
                         "Sum of matrices with different sizes lhs:[%zu][%zu], rhs:[%zu][%zu]",
                         amount_of_rows_, amount_of_columns_, rhs.amount_of_rows_, rhs.amount_of_columns_);
            throw std::runtime_error(error);
        }
        std::transform(data_.begin(), data_.end(), rhs.data_.begin(), data_.begin(),
                       [](double &lhs_value, double const &rhs_value) { return lhs_value += rhs_value; });
        return *this;
    }

    matrix operator+(matrix const &lhs, matrix const &rhs) {
        return matrix(lhs) += rhs;
    }

    matrix &matrix::operator-=(matrix const &rhs) {
        if (amount_of_rows_ != rhs.amount_of_rows_ || amount_of_columns_ != rhs.amount_of_columns_) {
            char error[256];
            std::sprintf(error,
                         "Sub of matrices with different sizes lhs:[%zu][%zu], rhs:[%zu][%zu]",
                         amount_of_rows_, amount_of_columns_, rhs.amount_of_rows_, rhs.amount_of_columns_);
            throw std::runtime_error(error);
        }
        std::transform(data_.begin(), data_.end(), rhs.data_.begin(), data_.begin(),
                       [](double &lhs_value, double const &rhs_value) { return lhs_value -= rhs_value; });
        return *this;
    }

    matrix operator-(matrix const &lhs, matrix const &rhs) {
        return matrix(lhs) -= rhs;
    }

    matrix &matrix::operator/=(double const &rhs) {
        if (rhs == 0) {
            throw std::runtime_error("Divizion by zero");
        }
        std::transform(data_.begin(), data_.end(), data_.begin(),
                       [rhs](double &lhs_value) { return lhs_value /= rhs; });
        return *this;
    }

    matrix operator/(matrix const &lhs, double const &rhs) {
        return matrix(lhs) /= rhs;
    }

    matrix &matrix::operator*=(double const &rhs) {
        if (rhs == 0) {
            throw std::runtime_error("Divizion by zero");
        }
        std::transform(data_.begin(), data_.end(), data_.begin(),
                       [rhs](double &lhs_value) { return lhs_value *= rhs; });
        return *this;
    }

    matrix operator*(matrix const &lhs, double const &rhs) {
        return matrix(lhs) *= rhs;
    }

    matrix operator*(double const &lhs, matrix const &rhs) {
        return matrix(rhs) *= lhs;
    }

    matrix &matrix::operator*=(matrix const &rhs) {
        if (amount_of_columns_ != rhs.amount_of_rows_) {
            char error[256];
            std::sprintf(error,
                         "Mul of matrices with bad sizes lhs:[%zu][%zu], rhs:[%zu][%zu]",
                         amount_of_rows_, amount_of_columns_, rhs.amount_of_rows_, rhs.amount_of_columns_);
            throw std::runtime_error(error);
        }
        matrix tmp(amount_of_rows_, rhs.amount_of_columns_);
        for (size_t i = 0; i < amount_of_rows_; i++) {
            for (size_t j = 0; j < rhs.amount_of_columns_; j++) {
                tmp(i, j) = operator()(i, 0) * rhs(0, j);
                for (size_t k = 1; k < amount_of_columns_; k++) {
                    tmp(i, j) += operator()(i, k) * rhs(k, j);
                }
            }
        }
        this->data_ = tmp.data_;
        amount_of_columns_ = rhs.amount_of_columns_;
        return *this;
    }

    matrix operator*(matrix const &lhs, matrix &rhs) {
        return matrix(lhs) *= rhs;
    }

    bool operator==(matrix const &lhs, matrix const &rhs) {
        if (lhs.amount_of_rows_ != rhs.amount_of_rows_ ||
            lhs.amount_of_columns_ != rhs.amount_of_columns_)
            return false;
        return std::equal(lhs.data_.begin(), lhs.data_.end(), rhs.data_.begin());
    }

    bool operator!=(matrix const &lhs, matrix const &rhs) {
        return !(lhs == rhs);
    }

    double matrix::determinant() const {
        if (amount_of_rows_ != amount_of_columns_) {
            char error[256];
            std::sprintf(error,
                         "Can't get determinant from non-square matrix [%zu][%zu]",
                         amount_of_rows_, amount_of_columns_);
            throw std::runtime_error(error);
        }

        std::vector<double> tmp(data_.begin(), data_.end());
        double res = 1.0;
        for (size_t i = 0; i < amount_of_rows_; i++) {
            size_t change = i;
            for (; change < amount_of_rows_; change++) {
                if (tmp[change * amount_of_columns_ + i] != 0) break;
            }
            if (change == amount_of_rows_) {
                return 0.0;
            }

            if (change != i) {
                res *= -1.0;
                std::swap_ranges(tmp.begin() + i * amount_of_columns_,
                                 tmp.begin() + (i + 1) * amount_of_columns_,
                                 tmp.begin() + change * amount_of_columns_);
            }
            for (size_t k = i + 1; k < amount_of_rows_; k++) {
                double mn = -1 * (tmp[k * amount_of_columns_ + i] / tmp[i * amount_of_columns_ + i]);
                std::transform(tmp.begin() + i * amount_of_columns_,
                               tmp.begin() + (i + 1) * amount_of_columns_,
                               tmp.begin() + k * amount_of_columns_,
                               tmp.begin() + k * amount_of_columns_,
                               [mn](double lhs, double rhs) { return rhs + mn * lhs; });
            }
        }

        for (size_t i = 0; i < amount_of_rows_; i++)
            res *= tmp[i * amount_of_columns_ + i];
        return res;
    }

    matrix &matrix::transpose() {
        if (amount_of_rows_ == 1 || amount_of_columns_ == 1) {
            std::swap(amount_of_rows_, amount_of_columns_);
            return *this;
        }
        if (amount_of_rows_ == amount_of_columns_) {
            for (size_t i = 0; i < amount_of_rows_; i++) {
                for (size_t j = 0; j < amount_of_columns_; j++) {
                    std::swap(data_[i * amount_of_columns_ + j], data_[i + j * amount_of_rows_]);
                }
            }
            return *this;
        }
        std::vector<bool> state(amount_of_columns_ * amount_of_rows_, false);
        size_t i = 0;
        while (i < data_.size()) {
            if (state[i]) {
                i++;
                continue;
            }
            state[i] = true;
            size_t save_state = i;
            double tmp = data_[i];
            do {
                i = (i % amount_of_columns_) * amount_of_rows_ + i / amount_of_columns_;
                std::swap(tmp, data_[i]);
                state[i] = true;
            } while (i != save_state);
        }

        std::swap(amount_of_rows_, amount_of_columns_);

        return *this;
    }

    matrix &matrix::inverse() {
        if (amount_of_rows_ != amount_of_columns_) {
            char error[256];
            std::sprintf(error,
                         "Can't inverse non-square matrix [%zu][%zu]",
                         amount_of_rows_, amount_of_columns_);
            throw std::runtime_error(error);
        }

        std::vector<double> tmp[2] = {std::vector<double>(amount_of_rows_ * amount_of_rows_),
                                      std::vector<double>(amount_of_rows_ * amount_of_rows_)};
        for (size_t i = 0; i < data_.size(); i++) {
            tmp[0][i] = data_[i];
        }
        for (size_t i = 0; i < amount_of_rows_; i++) {
            for (size_t j = 0; j < amount_of_columns_; j++) {
                tmp[1][i * amount_of_columns_ + j] = (i == j) ? 1.0 : 0.0;
            }
        }

        for (size_t i = 0; i < amount_of_rows_; i++) {
            size_t change = i;
            for (; change < amount_of_rows_; change++) {
                if (tmp[0][change * amount_of_columns_ + i] != 0) break;
            }
            if (change == amount_of_rows_) {
                throw std::runtime_error("Can't inverse matrix matrix with determinant = 0");
            }

            if (change != i) {
                for (auto &t : tmp) {
                    std::swap_ranges(t.begin() + i * amount_of_columns_,
                                     t.begin() + (i + 1) * amount_of_columns_,
                                     t.begin() + change * amount_of_columns_);
                }
            }
            double mn = 1.0 / tmp[0][i * amount_of_columns_ + i];
            for (auto &t : tmp) {
                std::for_each(t.begin() + i * amount_of_columns_,
                              t.begin() + (i + 1) * amount_of_columns_,
                              [mn](auto &lhs) { lhs *= mn; });
            }

            for (size_t k = i + 1; k < amount_of_rows_; k++) {
                mn = -1.0 * tmp[0][k * amount_of_columns_ + i];
                for (auto &t : tmp) {
                    std::transform(t.begin() + i * amount_of_columns_,
                                   t.begin() + (i + 1) * amount_of_columns_,
                                   t.begin() + k * amount_of_columns_,
                                   t.begin() + k * amount_of_columns_,
                                   [mn](double lhs, double rhs) { return rhs + mn * lhs; });
                }
            }
        }
        for (size_t i = amount_of_rows_; i--;) {
            for (size_t j = i + 1; j < amount_of_rows_; j++) {
                double mn = -1.0 * tmp[0][i * amount_of_columns_ + j];
                std::transform(tmp[1].begin() + i * amount_of_columns_,
                               tmp[1].begin() + (i + 1) * amount_of_columns_,
                               tmp[1].begin() + j * amount_of_columns_,
                               tmp[1].begin() + i * amount_of_columns_,
                               [mn](double lhs, double rhs) { return lhs + mn * rhs; });
            }
        }
        this->data_ = std::vector<double>(tmp[1].begin(), tmp[1].end());
        return *this;
    }

    size_t matrix::height() const {
        return amount_of_rows_;
    }

    size_t matrix::length() const {
        return amount_of_columns_;
    }

    matrix transpose(matrix const &a) {
        return matrix(a).transpose();
    }

    matrix inverse(matrix const &a) {
        return matrix(a).inverse();
    }

    size_t tensor::binary_pow(size_t const &dim, size_t const &power) {
        if (power == 0) {
            return 1;
        }
        if (power == 1) {
            return dim;
        }
        size_t tmp = binary_pow(dim, power / 2);
        if (power % 2 == 0) {
            return tmp * tmp;
        }
        return tmp * tmp * dim;
    }

    double tensor::sum_el(std::vector<size_t> const &old, std::vector<size_t> const &real, matrix T_matrix, matrix S_matrix) {
        double res = operator()(old);
        for (size_t i = 0; i < amount_of_p_; i++)
            res *= S_matrix(real[i], old[i]);
        for (size_t i = 0; i < amount_of_q_; i++)
            res *= T_matrix(old[i + amount_of_p_], real[i + amount_of_p_]);
        return res;
    }

    bool tensor::change_round(std::vector<size_t> &round) {
        size_t i = 0;
        while ((i != round.size()) && (round[i] + 1 == dimension_)) {
            round[i++] = 0;
        }
        if (i < round.size()) {
            round[i]++;
            return true;
        }
        return false;
    }

    double const &tensor::operator()(std::vector<size_t> const &v) const {
        if (v.size() != amount_of_p_ + amount_of_q_) {
            throw std::runtime_error("Can't get element from bad-size vector");
        }
        size_t real_pos = 0;
        size_t dim_tmp = 1;
        for (size_t i = 0; i < v.size(); i += 2) {
            if (v[i] >= dimension_) {
                char error[256];
                std::sprintf(error,
                             "Can't get element:[...][%zu][...], of tensor with dimension %zu",
                             v[i], dimension_);
                throw std::runtime_error(error);
            }
            if (v.size() - i == 1) {
                real_pos += v[i] * dim_tmp;
                continue;
            }
            real_pos += v[i + 1] * dim_tmp;
            dim_tmp *= dimension_;
            real_pos += v[i] * dim_tmp;
            dim_tmp *= dimension_;
        }
        return data_[real_pos];
    }

    double &tensor::operator()(std::vector<size_t> const &v) {
        return const_cast<double &>(const_cast<const tensor *>(this)->operator()(v));
    }

    tensor &tensor::contraction(size_t p_num, size_t q_num) {
        if (p_num >= amount_of_p_ || q_num >= amount_of_q_) {
            char error[256];
            std::sprintf(error,
                         "Can't contract tensor (p:%zu,q:%zu) by (%zu, %zu)",
                         amount_of_p_, amount_of_q_, p_num, q_num);
            throw std::runtime_error(error);
        }
        tensor result(dimension_, amount_of_p_ - 1, amount_of_q_ - 1);
        std::vector<size_t> round(amount_of_p_ + amount_of_q_ - 2, 0);
        std::vector<size_t> round_real(amount_of_p_ + amount_of_q_);
        bool no_next = false;
        do {
            for (size_t i = 0, j = 0; i < round_real.size(); i++) {
                if (i == p_num || i == amount_of_p_ + q_num) continue;
                round_real[i] = round[j++];
            }
            round_real[p_num] = 0;
            round_real[amount_of_p_ + q_num] = 0;
            result(round) = operator()(round_real);
            for (size_t i = 1; i < dimension_; i++) {
                round_real[p_num] = i;
                round_real[amount_of_p_ + q_num] = i;
                result(round) += operator()(round_real);
            }
        } while (change_round(round));
        *this = result;
        return *this;
    }

    tensor tensor::change_basis(matrix const &T_matrix) {
        if (T_matrix.height() != dimension_ || T_matrix.length() != dimension_) {
            char error[256];
            std::sprintf(error,
                         "Can't use double matrix[%zu][%zu] for tensor with dimension %zu",
                         T_matrix.height(), T_matrix.length(), dimension_);
            throw std::runtime_error(error);
        }
        matrix S_matrix = transpose(T_matrix);

        std::vector<size_t> round_real(amount_of_p_ + amount_of_q_, 0);
        std::vector<size_t> round_old(amount_of_p_ + amount_of_q_);

        bool no_next = false;
        bool no_next_old;
        size_t i = 0;
        tensor result(dimension_, amount_of_p_, amount_of_q_);
        do {
            round_old.assign(round_old.size(), 0);
            result(round_real) = sum_el(round_old, round_real, T_matrix, S_matrix);
            while (change_round(round_old)) {
                result(round_real) += sum_el(round_old, round_real, T_matrix, S_matrix);
            }
        } while (change_round(round_real));
        return result;
    }

    tensor contraction(tensor const &a, size_t p_num, size_t q_num) {
        return tensor(a).contraction(p_num, q_num);
    }
}