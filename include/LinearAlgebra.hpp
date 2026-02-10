#pragma once
#include <vector>
#include <cassert>

namespace LA {
    struct Matrix;

    auto HouseholderQR(Matrix& A, std::vector<double>& tau) -> void;

    auto computeHouseholder(Matrix& A, size_t k, std::vector<double>& tau) -> void;

    auto applyHouseholdToCol(Matrix& A, size_t k, size_t j, double tau) -> void;

    auto inPlaceQT(const Matrix& A, const std::vector<double>& tau, Matrix& Y) -> void;

    auto backSub(const Matrix& A, Matrix& B) -> void;
}
