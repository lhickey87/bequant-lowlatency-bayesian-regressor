#include "Regression.hpp"


auto OLS::FitInPlace(LA::Matrix& X, LA::Column y, LA::Column beta) -> void
{
    LA::QRSolve(X, y, beta);
}

auto OLS::Predict(const LA::Matrix &X, LA::ConstColumn beta) -> std::vector<double>
{
    return X*beta;
}
