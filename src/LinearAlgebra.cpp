#include "LinearAlgebra.hpp"
#include "Matrix.hpp"

namespace LA {

    auto computeHouseholder(Matrix& A, size_t k, std::vector<double>& tau) -> void
    {
        const double x0 = A(k, k);
        const auto A_K = A.SubColumn(k);

        const auto tailSS = A_K*A_K;

        if (tailSS == 0.0) { return;}

        const double normx = std::sqrt(x0 * x0 + tailSS);
        const double diag = (x0 >= 0.0) ? -normx : normx;

        tau[k] = (diag - x0) / diag;

        const double denom = (x0 - diag);
        assert(denom != 0.0);

        for (size_t i = k + 1; i < A.rows(); ++i) {
            A(i, k) /= denom;
        }
        A(k, k) = diag;
    }


    void HouseholderQR(Matrix& A, std::vector<double>& tau)
    {
        const size_t rows = A.rows();
        const size_t cols = A.cols();
        assert(rows >= cols);
        tau.assign(cols, 0.0);

        for (size_t k = 0; k < cols; ++k) {
            computeHouseholder(A, k, tau);

            const double tk = tau[k];
            if (tk == 0.0) continue;

            for (size_t j = k + 1; j < cols; ++j) {
                // A = A - tau * v * (v^T * A)
                applyHouseholdToCol(A, k, j, tk);
            }
        }
    }

    void backSub(const Matrix& R, Matrix& B)
    {
        const size_t n = R.cols();
        assert(R.rows() >= n);
        assert(B.rows() == n);

        const size_t rhsCount = B.cols();

        for (size_t i = n; i-- > 0;) {
            const double rii = R(i, i);
            assert(rii != 0.0);

            for (size_t col = 0; col < rhsCount; ++col) {
                double sum = B(i, col);
                for (size_t j = i + 1; j < n; ++j) {
                    sum -= R(i, j) * B(j, col);
                }
                B(i, col) = sum / rii;
            }
        }
    }

    void applyHouseholdToCol(Matrix& A, size_t k, size_t j, double tau)
    {
        auto A_k = A.SubColumn(k);
        auto A_j = A.SubColumn(j);

        double dotProd = (A(k, j)+A_k*A_j)*tau;

        A(k, j) -= dotProd; // v0=1

        A_j.addScaled(A_k, -dotProd);
    }


    void inPlaceQT(Matrix& A, const std::vector<double>& tau, Matrix& y){
        const size_t m = A.rows();
        const size_t n = A.cols();

        assert(y.rows() == m && y.cols() == 1);
        assert(tau.size() == n);

        for (size_t k = 0; k < n; ++k) {
            const double tk = tau[k];

            if (tk == 0.0) continue;

            auto A_K = A.SubColumn(k);
            auto Y = y.getColumn(0);

            double dot = (y(k, 0) + Y*A_K)*tk;              // v0 = 1
            y(k, 0) -= dot;

            //Y = Y - A_K*dot
            Y.addScaled(A_K,-dot);
        }
    }

}
