#pragma once
#include <span>
#include <omp.h>
#include <cassert>
#include "LinearAlgebra.hpp"
#include "ColumnView.hpp"

namespace LA {

    struct Matrix {

        Matrix() = delete;
        explicit Matrix(size_t rows, size_t cols);
        explicit Matrix(size_t rows, size_t cols, std::span<const double> data);

        inline double& operator()(size_t row, size_t col) noexcept {
            return data_[col*rows_+cols_];
        }

        inline const double& operator()(size_t row, size_t col) const noexcept {
            return data_[col*rows_+cols_];
        }

        auto operator*(const Matrix& rhs) const -> Matrix;

        auto inline rows() const noexcept -> size_t { return rows_;}
        auto inline cols() const noexcept -> size_t {return cols_;}

        auto SubColumn(size_t column) const -> ColumnView;
        auto getColumn(size_t column) const -> ColumnView;

        auto transpose() const -> Matrix;

        auto Solve(const Matrix& mat) const -> Matrix;
        auto QRSolve(const Matrix& y) const -> Matrix;
    private:

        std::vector<double> data_;
        size_t rows_{};
        size_t cols_{};
    };

}
