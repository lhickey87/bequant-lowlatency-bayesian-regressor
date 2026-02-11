#pragma once
#include "Matrix.hpp"
#include "ColumnView.hpp"
#include <iostream>

namespace OLS {
    auto FitInPlace(LA::Matrix& X, LA::Column Y, LA::Column beta) -> void;
    auto Predict(const LA::Matrix& X, LA::ConstColumn beta) -> std::vector<double>;
};
