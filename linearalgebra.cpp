#include "linearalgebra.h"

namespace la {
double const &matrix::operator()(size_t const &row, size_t const &column) const {
    if (row >= amount_of_rows_ || column >= amount_of_columns_) {
        char error[256];
        std::sprintf(error, "Adress to [%zu][%zu], size of matrix is [%zu][%zu]", row, column, amount_of_rows_,
                     amount_of_columns_);
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
        std::sprintf(error, "Sum of matrices with different sizes lhs:[%zu][%zu], rhs:[%zu][%zu]", amount_of_rows_,
                     amount_of_columns_, rhs.amount_of_rows_, rhs.amount_of_columns_);
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
        std::sprintf(error, "Sub of matrices with different sizes lhs:[%zu][%zu], rhs:[%zu][%zu]", amount_of_rows_,
                     amount_of_columns_, rhs.amount_of_rows_, rhs.amount_of_columns_);
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
    std::transform(data_.begin(), data_.end(), data_.begin(), [rhs](double &lhs_value) { return lhs_value /= rhs; });
    return *this;
}

matrix operator/(matrix const &lhs, double const &rhs) {
    return matrix(lhs) /= rhs;
}

matrix &matrix::operator*=(double const &rhs) {
    if (rhs == 0) {
        throw std::runtime_error("Divizion by zero");
    }
    std::transform(data_.begin(), data_.end(), data_.begin(), [rhs](double &lhs_value) { return lhs_value *= rhs; });
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
        std::sprintf(error, "Mul of matrices with bad sizes lhs:[%zu][%zu], rhs:[%zu][%zu]", amount_of_rows_,
                     amount_of_columns_, rhs.amount_of_rows_, rhs.amount_of_columns_);
        throw std::runtime_error(error);
    }
    matrix tmp(amount_of_rows_, rhs.amount_of_columns_);
    for (size_t i = 0; i < amount_of_rows_; ++i) {
        for (size_t j = 0; j < rhs.amount_of_columns_; ++j) {
            tmp(i, j) = operator()(i, 0) * rhs(0, j);
            for (size_t k = 1; k < amount_of_columns_; ++k) {
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
    if (lhs.amount_of_rows_ != rhs.amount_of_rows_ || lhs.amount_of_columns_ != rhs.amount_of_columns_)
        return false;
    return std::equal(lhs.data_.begin(), lhs.data_.end(), rhs.data_.begin());
}

bool operator!=(matrix const &lhs, matrix const &rhs) {
    return !(lhs == rhs);
}

double matrix::determinant() const {
    if (amount_of_rows_ != amount_of_columns_) {
        char error[256];
        std::sprintf(error, "Can't get determinant from non-square matrix [%zu][%zu]", amount_of_rows_, amount_of_columns_);
        throw std::runtime_error(error);
    }

    std::vector<double> tmp(data_.begin(), data_.end());
    double res = 1.0;
    for (size_t i = 0; i < amount_of_rows_; ++i) {
        size_t change = i;
        for (; change < amount_of_rows_; ++change) {
            if (tmp[change * amount_of_columns_ + i] != 0)
                break;
        }
        if (change == amount_of_rows_) {
            return 0.0;
        }

        if (change != i) {
            res *= -1.0;
            std::swap_ranges(tmp.begin() + i * amount_of_columns_, tmp.begin() + (i + 1) * amount_of_columns_,
                             tmp.begin() + change * amount_of_columns_);
        }
        for (size_t k = i + 1; k < amount_of_rows_; ++k) {
            double mn = -1 * (tmp[k * amount_of_columns_ + i] / tmp[i * amount_of_columns_ + i]);
            std::transform(tmp.begin() + i * amount_of_columns_, tmp.begin() + (i + 1) * amount_of_columns_,
                           tmp.begin() + k * amount_of_columns_, tmp.begin() + k * amount_of_columns_,
                           [mn](double lhs, double rhs) { return rhs + mn * lhs; });
        }
    }

    for (size_t i = 0; i < amount_of_rows_; ++i)
        res *= tmp[i * amount_of_columns_ + i];
    return res;
}

matrix &matrix::transpose() {
    if (amount_of_rows_ == 1 || amount_of_columns_ == 1) {
        std::swap(amount_of_rows_, amount_of_columns_);
        return *this;
    }
    if (amount_of_rows_ == amount_of_columns_) {
        for (size_t i = 0; i < amount_of_rows_ / 2; ++i) {
            for (size_t j = 0; j < amount_of_columns_ / 2; ++j) {
                std::swap(data_[i * amount_of_columns_ + j], data_[i + j * amount_of_rows_]);
            }
        }
        return *this;
    }
    std::vector<bool> state(amount_of_columns_ * amount_of_rows_, false);
    size_t i = 0;
    while (i < data_.size()) {
        if (state[i]) {
            ++i;
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
        std::sprintf(error, "Can't inverse non-square matrix [%zu][%zu]", amount_of_rows_, amount_of_columns_);
        throw std::runtime_error(error);
    }

    std::vector<double> tmp[2] = {std::vector<double>(amount_of_rows_ * amount_of_rows_),
                                  std::vector<double>(amount_of_rows_ * amount_of_rows_)};
    for (size_t i = 0; i < data_.size(); ++i) {
        tmp[0][i] = data_[i];
    }
    for (size_t i = 0; i < amount_of_rows_; ++i) {
        for (size_t j = 0; j < amount_of_columns_; ++j) {
            tmp[1][i * amount_of_columns_ + j] = (i == j) ? 1.0 : 0.0;
        }
    }

    for (size_t i = 0; i < amount_of_rows_; ++i) {
        size_t change = i;
        for (; change < amount_of_rows_; ++change) {
            if (tmp[0][change * amount_of_columns_ + i] != 0)
                break;
        }
        if (change == amount_of_rows_) {
            throw std::runtime_error("Can't inverse matrix with determinant = 0");
        }

        if (change != i) {
            for (auto &t : tmp) {
                std::swap_ranges(t.begin() + i * amount_of_columns_, t.begin() + (i + 1) * amount_of_columns_,
                                 t.begin() + change * amount_of_columns_);
            }
        }
        double mn = 1.0 / tmp[0][i * amount_of_columns_ + i];
        for (auto &t : tmp) {
            std::for_each(t.begin() + i * amount_of_columns_, t.begin() + (i + 1) * amount_of_columns_,
                          [mn](auto &lhs) { lhs *= mn; });
        }

        for (size_t k = i + 1; k < amount_of_rows_; ++k) {
            mn = -1.0 * tmp[0][k * amount_of_columns_ + i];
            for (auto &t : tmp) {
                std::transform(t.begin() + i * amount_of_columns_, t.begin() + (i + 1) * amount_of_columns_,
                               t.begin() + k * amount_of_columns_, t.begin() + k * amount_of_columns_,
                               [mn](double lhs, double rhs) { return rhs + mn * lhs; });
            }
        }
    }
    for (size_t i = amount_of_rows_; --i;) {
        for (size_t j = i + 1; j < amount_of_rows_; ++j) {
            double mn = -1.0 * tmp[0][i * amount_of_columns_ + j];
            std::transform(tmp[1].begin() + i * amount_of_columns_, tmp[1].begin() + (i + 1) * amount_of_columns_,
                           tmp[1].begin() + j * amount_of_columns_, tmp[1].begin() + i * amount_of_columns_,
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

std::ostream &operator<<(std::ostream &cout, matrix const &a) {
    for (size_t i = 0; i < a.amount_of_columns_; ++i) {
        for (size_t j = 0; j < a.amount_of_rows_; ++j) {
            cout << a(i, j) << ' ';
        }
        cout << '\n';
    }
    return cout;
}

std::istream &operator>>(std::istream &cin, matrix &a) {
    size_t h, l;
    cin >> h >> l;
    a = matrix(h, l);
    for (size_t i = 0; i < h; ++i) {
        for (size_t j = 0; j < l; ++j) {
            cin >> a(i, j);
        }
    }
    return cin;
}

matrix transpose(matrix const &a) {
    return matrix(a).transpose();
}

matrix inverse(matrix const &a) {
    return matrix(a).inverse();
}

matrix E(size_t n) {
    matrix result(n, n);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            result(i, j) = (i == j) ? 1.0 : 0.0;
        }
    }
    return result;
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
    for (size_t i = 0; i < amount_of_p_; ++i)
        res *= S_matrix(real[i], old[i]);
    for (size_t i = 0; i < amount_of_q_; ++i)
        res *= T_matrix(old[i + amount_of_p_], real[i + amount_of_p_]);
    return res;
}

bool tensor::change_round(std::vector<size_t> &round) const {
    size_t i = 0;
    size_t real_i = (round.size() > 1) ? 1 : 0;
    while ((i != round.size()) && (round[real_i] + 1 == dimension_)) {
        round[real_i] = 0;
        ++i;
        real_i = (round.size() % 2 == 0 || round.size() - 1 != i) ? i ^ 1 : i;
    }
    if (i < round.size()) {
        ++round[real_i];
        return true;
    }
    return false;
}

bool next_permutation_mn(std::vector<size_t> &permutation, double &mn) {
    size_t i = permutation.size() - 1;
    while (i > 0 && permutation[i] < permutation[i - 1]) {
        --i;
    }

    if (i == 0) {
        return false;
    }
    size_t pos = i;
    for (size_t j = i + 1; j < permutation.size(); ++j) {
        if (permutation[j] > permutation[i - 1] && permutation[j] < permutation[pos]) {
            pos = j;
        }
    }
    std::swap(permutation[pos], permutation[i - 1]);
    mn *= -1.0;

    bool is_sorted = false;
    while (!is_sorted) {
        is_sorted = true;
        for (size_t j = i + 1; j < permutation.size(); ++j) {
            if (permutation[j] < permutation[j - 1]) {
                is_sorted = false;
                mn *= -1.0;
                std::swap(permutation[j], permutation[j - 1]);
            }
        }
    }
    return true;
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
            std::sprintf(error, "Can't get element:[...][%zu][...], of tensor with dimension %zu", v[i], dimension_);
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

double const &tensor::operator()() const {
    if (amount_of_p_ + amount_of_q_ > 0) {
        throw std::runtime_error("Use operator (vector) instead of operator () for non-scalar tensor");
    }
    return data_.front();
}

double &tensor::operator()() {
    return const_cast<double &>(const_cast<const tensor *>(this)->operator()());
}

tensor &tensor::contraction(size_t p_num, size_t q_num) {
    if (p_num >= amount_of_p_ || q_num >= amount_of_q_) {
        char error[256];
        std::sprintf(error, "Can't contract tensor (p:%zu,q:%zu) by (%zu, %zu)", amount_of_p_, amount_of_q_, p_num, q_num);
        throw std::runtime_error(error);
    }
    tensor result(dimension_, amount_of_p_ - 1, amount_of_q_ - 1);
    std::vector<size_t> round(amount_of_p_ + amount_of_q_ - 2, 0);
    std::vector<size_t> round_real(amount_of_p_ + amount_of_q_);
    do {
        for (size_t i = 0, j = 0; i < round_real.size(); ++i) {
            if (i == p_num || i == amount_of_p_ + q_num)
                continue;
            round_real[i] = round[j++];
        }
        round_real[p_num] = 0;
        round_real[amount_of_p_ + q_num] = 0;
        result(round) = operator()(round_real);
        for (size_t i = 1; i < dimension_; ++i) {
            round_real[p_num] = i;
            round_real[amount_of_p_ + q_num] = i;
            result(round) += operator()(round_real);
        }
    } while (change_round(round));
    *this = result;
    return *this;
}

tensor &tensor::change_basis(matrix const &T_matrix) {
    if (T_matrix.height() != dimension_ || T_matrix.length() != dimension_) {
        char error[256];
        std::sprintf(error, "Can't use double matrix[%zu][%zu] for tensor with dimension %zu", T_matrix.height(),
                     T_matrix.length(), dimension_);
        throw std::runtime_error(error);
    }
    matrix S_matrix = inverse(T_matrix);

    std::vector<size_t> round_real(amount_of_p_ + amount_of_q_, 0);
    std::vector<size_t> round_old(amount_of_p_ + amount_of_q_);

    size_t i = 0;
    tensor result(dimension_, amount_of_p_, amount_of_q_);
    do {
        round_old.assign(round_old.size(), 0);
        result(round_real) = sum_el(round_old, round_real, T_matrix, S_matrix);
        while (change_round(round_old)) {
            result(round_real) += sum_el(round_old, round_real, T_matrix, S_matrix);
        }
    } while (change_round(round_real));
    *this = result;
    return *this;
}

tensor &tensor::sym(std::vector<bool> const &v) {
    if (v.empty() && amount_of_p_ + amount_of_q_ == 0) {
        return *this;
    }
    if (v.size() != amount_of_p_ + amount_of_q_) {
        throw std::runtime_error("Can't get element from bad-size vector");
    }
    size_t cnt = 0;
    for (size_t i = 0; i < v.size(); ++i) {
        if (v[i]) {
            ++cnt;
        }
    }
    if (cnt < 2) {
        return *this;
    }

    tensor tmp(dimension_, amount_of_p_, amount_of_q_);
    tensor state(dimension_, amount_of_p_, amount_of_q_);

    for (double &i : state.data_) {
        i = -1;
    }

    std::vector<size_t> round(amount_of_p_ + amount_of_q_);
    std::vector<size_t> state_round(amount_of_p_ + amount_of_q_, 0);
    do {
        if (state(state_round) > 0)
            continue;
        double sum = 0;
        size_t amount = 0;
        std::vector<size_t> permutation(cnt);
        for (size_t i(0), j(0); i < v.size(); ++i) {
            if (!v[i])
                continue;
            permutation[j++] = i;
        }
        do {
            for (size_t i(0), j(0); i < round.size(); ++i) {
                if (v[i]) {
                    round[i] = state_round[permutation[j++]];
                } else {
                    round[i] = state_round[i];
                }
            }
            ++amount;
            sum += operator()(round);
        } while (std::next_permutation(permutation.begin(), permutation.end()));
        sum /= amount;
        for (size_t i(0), j(0); i < v.size(); ++i) {
            if (!v[i])
                continue;
            permutation[j++] = i;
        }
        do {
            for (size_t i(0), j(0); i < round.size(); ++i) {
                if (v[i]) {
                    round[i] = state_round[permutation[j++]];
                } else {
                    round[i] = state_round[i];
                }
            }
            tmp(round) = sum;
            state(round) = 1;
        } while (std::next_permutation(permutation.begin(), permutation.end()));
    } while (change_round(state_round));
    *this = tmp;
    return *this;
}

tensor &tensor::asym(std::vector<bool> const &v) {
    if (v.empty() && amount_of_p_ + amount_of_q_ == 0) {
        return *this;
    }
    if (v.size() != amount_of_p_ + amount_of_q_) {
        throw std::runtime_error("Can't get element from bad-size vector");
    }
    size_t cnt = 0;
    for (size_t i = 0; i < v.size(); ++i) {
        if (v[i]) {
            ++cnt;
        }
    }
    if (cnt < 2) {
        return *this;
    }

    tensor tmp(dimension_, amount_of_p_, amount_of_q_);
    tensor state(dimension_, amount_of_p_, amount_of_q_);

    for (double &i : state.data_) {
        i = -1;
    }

    std::vector<size_t> round(amount_of_p_ + amount_of_q_);
    std::vector<size_t> state_round(amount_of_p_ + amount_of_q_, 0);
    do {
        if (state(state_round) > 0)
            continue;
        double mn = 1.0;
        double sum = 0;
        size_t amount = 0;
        std::vector<size_t> permutation(cnt);
        for (size_t i(0), j(0); i < v.size(); ++i) {
            if (!v[i])
                continue;
            permutation[j++] = i;
        }
        do {
            for (size_t i(0), j(0); i < round.size(); ++i) {
                if (v[i]) {
                    round[i] = state_round[permutation[j++]];
                } else {
                    round[i] = state_round[i];
                }
            }
            ++amount;
            sum += mn * operator()(round);
        } while (next_permutation_mn(permutation, mn));
        sum /= amount;
        for (size_t i(0), j(0); i < v.size(); ++i) {
            if (!v[i])
                continue;
            permutation[j++] = i;
        }
        mn = 1.0;
        do {
            for (size_t i(0), j(0); i < round.size(); ++i) {
                if (v[i]) {
                    round[i] = state_round[permutation[j++]];
                } else {
                    round[i] = state_round[i];
                }
            }
            tmp(round) = sum * mn;
            state(round) = 1;
        } while (next_permutation_mn(permutation, mn));
    } while (change_round(state_round));
    *this = tmp;
    return *this;
}

std::ostream &operator<<(std::ostream &cout, tensor const &a) {
    if (a.amount_of_p_ + a.amount_of_q_ == 0) {
        cout << a();
        return cout;
    }
    std::vector<size_t> v(a.amount_of_p_ + a.amount_of_q_, 0);
    std::vector<size_t> prev(a.amount_of_p_ + a.amount_of_q_, 0);
    bool first = true;
    do {
        size_t max = 0;
        for (size_t i = 0; i < v.size(); ++i) {
            if (v[i] != prev[i]) {
                max = std::max(max, (v.size() % 2 == 0 || v.size() - 1 != i) ? i ^ 1 : i);
            }
        }
        if (first) {
            first = false;
        } else if (max == 0) {
            cout << ' ';
        } else if (max == 1) {
            cout << '\n';
        } else {
            cout << '\n';
            for (size_t j = 0; j < max - 1; ++j)
                cout << "_____________\n";
        }
        cout << a(v);
        prev = v;
    } while (a.change_round(v));
    return cout;
}

std::istream &operator>>(std::istream &cin, tensor &a) {
    size_t d, p, q;
    cin >> d >> p >> q;
    a = tensor(d, p, q);
    if (p + q == 0) {
        cin >> a();
        return cin;
    }
    std::vector<size_t> v(p + q, 0);
    do {
        cin >> a(v);
    } while (a.change_round(v));
    return cin;
}

tensor contraction(tensor const &a, size_t p_num, size_t q_num) {
    return tensor(a).contraction(p_num, q_num);
}

tensor change_basis(tensor const &a, matrix const &T_matrix) {
    return tensor(a).change_basis(T_matrix);
}

tensor sym(tensor const &a, std::vector<bool> const &v) {
    return tensor(a).sym(v);
}

tensor asym(tensor const &a, std::vector<bool> const &v) {
    return tensor(a).sym(v);
}
}  // namespace la
