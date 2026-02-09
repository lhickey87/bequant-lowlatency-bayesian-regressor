#pragma once
#include <span>
#include <omp.h>
#include <cassert>
#include "LinearAlgebra.hpp"

struct Matrix {

    Matrix() = delete;
    explicit Matrix(size_t rows, size_t cols);
    explicit Matrix(std::span<double> data);

    Matrix(const Matrix&) = default;
    Matrix(Matrix&&) = default;

    inline double& operator()(size_t row, size_t col) noexcept;

    inline const double& operator()(size_t row, size_t col) const noexcept;

    Matrix operator*(Matrix& rhs);

    Matrix transpose() const;

    Matrix Solve(const Matrix& mat) const;

    Matrix QRSolve(const Matrix& y) const;
private:
    void HouseholderQR(Matrix& A, std::vector<double>& tau);

    inline void applyHouseholdToCol(Matrix& A, size_t k, size_t j, double tau);
    inline double tailSumSquares(Matrix& A, size_t startInd);

    void inPlaceQT(const Matrix& A, const std::vector<double>& tau, Matrix& Y);

    void backSub(const Matrix& A, Matrix& B);

    std::vector<double> data_;
    size_t rows_{};
    size_t cols_{};
};
