#include "Matrix.hpp"

namespace LA{

    Matrix::Matrix(size_t row, size_t col) : data_(row * col, 0.0), rows_(row), cols_(col) {}

    Matrix::Matrix(size_t rows, size_t cols, std::span<const double> data)
        : data_(data.begin(), data.end()), rows_(rows), cols_(cols) {
        assert(rows_ * cols_ == data_.size());
    }

    auto Matrix::transpose() const -> Matrix {
        Matrix result(cols_,rows_);
        //so now we should do most cache friendly possible way via row major order
        for (size_t i = 0; i < rows_; ++i){
            for (size_t j = 0; j < cols_; ++j){
                result(j,i) = (*this)(i,j);
            }
        }
        return result;
    }

    auto Matrix::SubColumn(size_t column) -> Column
    {
        return ColumnView(&data_[column*rows_+(column+1)], column*rows_);
    }

    auto Matrix::SubColumn(size_t column) const -> ConstColumn
    {
        return ConstColumn(&data_[column*rows_+(column+1)], column*rows_);
    }

    auto Matrix::getColumn(size_t column) -> Column
    {
        return ColumnView(&data_[column*rows_], column*rows_);
    }

    auto Matrix::getColumn(size_t column) const -> ConstColumn
    {
        return ConstColumn(&data_[column*rows_], column*rows_);
    }

    auto Matrix::operator*(const Matrix& rhs) const -> Matrix {
        assert(cols_ == rhs.rows_);

        Matrix result(rows_, rhs.cols_);
        #pragma omp parralel for if (rows_ > 128)
        for (size_t i = 0; i < rows_; ++i) {
            for (size_t k = 0; k < cols_; ++k) {
                auto aik = (*this)(i, k);          // reuse this value
                for (size_t j = 0; j < rhs.cols_; ++j) {
                    result(i, j) += aik * rhs(k, j);
                }
            }
        }
        return result;
    }

    auto operator*(const Matrix& X, ConstColumn col) ->
    {
        const size_t vecLength = col.length();
        assert(X.cols() == vecLength);

        Matrix result(vecLength,1);
        for (size_t i = 0; i < vecLength;++i)
        {
            result(i,0) =
        }
    }

}
