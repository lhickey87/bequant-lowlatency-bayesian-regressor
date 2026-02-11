#include "Regression.hpp"


auto OLS::FitInPlace(LA::Matrix& X, LA::Column y, LA::Column beta) -> void
{
    LA::QRSolve(X, y, beta);
}

auto OLS::Predict(const LA::Matrix &X, LA::ConstColumn beta) -> LA::Column
{
    //all we need to do is perform y = XB, X*B

}
