#include "ColumnView.hpp"


namespace LA{

    ColumnView::ColumnView(std::span<double> col) : data_(col){}

    ColumnView::ColumnView(double* data, size_t length) : data_(data,length){}

    auto ColumnView::operator*=(double scalar) -> void
    {
        for (size_t i = 0; i < size();++i)
        {
            data_[i] *= scalar;
        }
    }

    auto ColumnView::operator-=(ColumnView rhs) -> void
    {
        assert(size() == rhs.size());
        for (size_t i = 0; i < size(); ++i)
        {
            data_[i] -= rhs[i];
        }
    }

    auto operator*(ColumnView lhs, ColumnView rhs) noexcept -> double
    {
        double sum = 0.0;
        for (size_t i = 0; i < lhs.size();++i)
        {
            sum += lhs[i]*rhs[i];
        }
        return sum;
    }
}
