#pragma once
#include <vector>
#include <cassert>
#include "ColumnView.hpp"

namespace LA {
    struct Matrix;

    auto HouseholderQR(Matrix& X, std::vector<double>& tau) -> void;

    auto computeHouseholder(Matrix& X, size_t k, std::vector<double>& tau) -> void;

    auto applyHouseholdToCol(Matrix& X, size_t k, size_t j, double tau) -> void;

    auto inPlaceQT(const Matrix& QT, const std::vector<double>& tau, Column Y) -> void;

    auto backSub(const Matrix& R, ConstColumn yTop, Column result) -> void;

    auto QRSolve(Matrix& X, Column y, Column beta_out) -> void;
}
