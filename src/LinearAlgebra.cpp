#include "LinearAlgebra.hpp"
#include "Matrix.hpp"

namespace LA {

    auto computeHouseholder(Matrix& X, size_t k, std::vector<double>& tau) -> void
    {
        const double x0 = X(k, k);
        const auto X_k = X.SubColumn(k);

        const auto tailSS = X_k*X_k;

        if (tailSS == 0.0) { return; } // tau[k] already zero-initialized

        const double normx = std::sqrt(x0 * x0 + tailSS);
        const double diag = (x0 >= 0.0) ? -normx : normx;

        tau[k] = (diag - x0) / diag;

        const double denom = (x0 - diag);
        assert(denom != 0.0);

        for (size_t i = k + 1; i < X.rows(); ++i) {
            X(i, k) /= denom;
        }
        X(k, k) = diag;
    }


    void HouseholderQR(Matrix& X, std::vector<double>& tau)
    {
        const size_t rows = X.rows();
        const size_t cols = X.cols();

        assert(rows >= cols);
        assert(tau.size() == cols);

        for (size_t k = 0; k < cols; ++k) {
            computeHouseholder(X, k, tau);

            const double tk = tau[k];
            if (tk == 0.0) continue;

            for (size_t j = k + 1; j < cols; ++j) {
                // X = X - tau * v * (v^T * X)
                applyHouseholdToCol(X, k, j, tk);
            }
        }
    }

    void backSub(const Matrix& R, ConstColumn rhs_top, Column beta_out)
    {
        const size_t n = R.cols();
        assert(R.rows() >= n);
        assert(rhs_top.length() == n);
        assert(beta_out.length() == n);

        for (size_t i = 0; i < n; ++i) {
            beta_out[i] = rhs_top[i];
        }

        for (size_t i = n; i-- > 0;) {
            const double rii = R(i, i);
            assert(rii != 0.0);

            double sum = beta_out[i];
            for (size_t j = i + 1; j < n; ++j) {
                sum -= R(i, j) * beta_out[j];
            }
            beta_out[i] = sum / rii;
        }
    }

    void applyHouseholdToCol(Matrix& X, size_t k, size_t j, double tau)
    {
        auto X_k = X.SubColumn(k);
        auto X_j = X.SubColumn(j);

        double vTx = (X(k, j)+X_k*X_j)*tau;

        X(k, j) -= vTx; // v0=1

        X_j.addScaled(X_k, -vTx);
    }


    void inPlaceQT(Matrix& QT, const std::vector<double>& tau, Column y){
        const size_t m = QT.rows();
        const size_t n = QT.cols();

        assert(y.length() == m);
        assert(tau.size() == n);

        for (size_t k = 0; k < n; ++k) {
            const double tk = tau[k];

            if (tk == 0.0) continue;
            auto X_k = QT.SubColumn(k);

            double vTy = (y[k] + y*X_k)*tk;              // v0 = 1
            y[k] -= vTy;

            // y = y - v * (v^T y) * tau
            y.addScaled(X_k,-vTy);
        }
    }

    auto QRSolve(Matrix& X, Column y, Column beta) -> void
    {
        const auto n = X.cols();

        assert(X.rows() == y.length());
        assert(beta.length() == n);

        std::vector<double> tau(n,0.0);
        HouseholderQR(X, tau);   // X -> R (upper triangle)
        inPlaceQT(X, tau, y);    // y -> Q^T y

        // Solve: R * beta = (Q^T y)_top_n
        auto y_top = ConstColumn(y.data(), n);

        backSub(X, y_top, beta);
    }
}
