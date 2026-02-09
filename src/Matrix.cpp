#include "Matrix.hpp"

Matrix::Matrix(size_t row, size_t col) : data_(row*col,0){}

Matrix::Matrix(std::span<double> data){
    for (size_t i = 0; i < data.size(); ++i){
        data_[i] = data[i];
    }
}

Matrix Matrix::transpose() const {
    Matrix result(cols_,rows_);
    //so now we should do most cache friendly possible way via row major order
    for (size_t i = 0; i < rows_; ++i){
        for (size_t j = 0; i < cols_; ++j){
            result(j,i) = (*this)(i,j);
        }
    }
    return result;
}

// thus we would do row*cols_+ col
inline double& Matrix::operator()(size_t row, size_t col) noexcept{
    return data_[row*cols_+col];
}

inline const double& Matrix::operator()(size_t row, size_t col) const noexcept{
    return data_[row*cols_+col];
}

void Matrix::CholeskySolve(Matrix& mat) noexcept {
    assert(rows_ == cols_);
    for (size_t i = 0; i < rows_; ++i){

    }
}

inline void Matrix::applyHouseholdToCol(Matrix& A, size_t k, size_t j, double tau)
{
    double dotProd = A(k, j); // v0 = 1
    for (size_t i = k + 1; i < A.rows_; ++i) {
        dotProd += A(i, k) * A(i, j);
    }
    dotProd *= tau;

    A(k, j) -= dotProd; // v0=1
    for (size_t i = k + 1; i < A.rows_; ++i) {
        A(i, j) -= A(i, k) * dotProd;
    }
}

inline double Matrix::tailSumSquares(Matrix& A, size_t c)
{
    const size_t m = A.rows_;
    double sigma = 0.0;

    for (size_t i = c + 1; i < m; ++i) {
        const double xi = A(i, c);
        sigma += xi * xi;
    }
    return sigma;
}

void Matrix::HouseholderQR(Matrix& A, std::vector<double>& tau)
{
    const size_t m = A.rows_;
    const size_t n = A.cols_;
    assert(m >= n);
    tau.assign(n, 0.0);

    for (size_t k = 0; k < n; ++k) {
        const double x0 = A(k, k);

        const double tailSS = tailSumSquares(A,k);
        if (tailSS == 0.0) {
            tau[k] = 0.0;
            continue;
        }

        const double normx = std::sqrt(x0 * x0 + tailSS);
        const double diag = (x0 >= 0.0) ? -normx : normx;

        const double tk = (diag - x0) / diag;
        tau[k] = tk;
        const double denom = (x0 - diag);
        assert(denom != 0.0);

        for (size_t i = k + 1; i < m; ++i) {
            A(i, k) /= denom;            // store v_tail in-place
        }

        A(k, k) = diag;

        for (size_t j = k + 1; j < n; ++j) {
            //A = A - tau*v(v^T*A)
            applyHouseholdToCol(A, k, j, tk);
        }
    }
}

void Matrix::backSub(const Matrix& A, Matrix& B) {
    const size_t n = A.cols_;
    assert(A.rows_ >= n);
    assert(B.rows_ == n);

    const size_t kRHS = B.cols_;

    for (size_t ii = n; ii-- > 0;) {
        const double rii = A(ii, ii);
        assert(rii != 0.0);

        for (size_t col = 0; col < kRHS; ++col) {
            double sum = B(ii, col);
            for (size_t j = ii + 1; j < n; ++j) {
                sum -= A(ii, j) * B(j, col);
            }
            B(ii, col) = sum / rii;
        }
    }
}

void Matrix::inPlaceQT(const Matrix& A, const std::vector<double>& tau, Matrix& y){
    const size_t m = A.rows_;
    const size_t n = A.cols_;
    assert(y.rows_ == m && y.cols_ == 1);
    assert(tau.size() == n);

    for (size_t k = 0; k < n; ++k) {
        const double tk = tau[k];
        if (tk == 0.0) continue;

        double dot = y(k, 0);              // v0 = 1
        for (size_t i = k + 1; i < m; ++i) {
            dot += A(i, k) * y(i, 0);      // A(i,k) stores v_tail
        }
        dot *= tk;

        // y(k:m-1) -= v * dot
        // we are able to assume that v0 = 1 as we will perform normalization
        y(k, 0) -= dot;                    // v0 = 1
        for (size_t i = k + 1; i < m; ++i) {
            y(i, 0) -= A(i, k) * dot;
        }
    }
}

Matrix Matrix::operator*(Matrix& rhs){
    assert(cols_ == rhs.rows_);

    Matrix result(rows_, rhs.cols_);
    #pragma omp parralel for if (rows_ > 128)
    for (int i = 0; i < rows_; ++i) {
        for (int k = 0; k < cols_; ++k) {
            auto aik = (*this)(i, k);          // reuse this value
            for (int j = 0; j < rhs.cols_; ++j) {
                result(i, j) += aik * rhs(k, j);
            }
        }
    }
    return result;
}
