#pragma once
#include <span>
#include <numeric>
#include <cassert>

namespace LA {
    struct ColumnView {

        ColumnView(std::span<double> col);
        explicit ColumnView(const double* data, size_t length);
        explicit ColumnView(double* data, size_t length);

        auto inline operator[](size_t index) -> double&
        {
            return data_[index];
        }

        auto operator-=(ColumnView rhs) -> void;
        auto operator*=(double scalar) -> void;

        auto inline begin() {return data_.begin();}
        auto inline end() {return data_.end();}
        auto inline size() const -> size_t{ return data_.size();}

        auto inline addScaled(ColumnView x, double a) noexcept
        {
            assert(size() == x.size());
            for (size_t i = 0; i < size(); ++i) {
                data_[i] += a * x[i];
            }
        }

        auto inline SubColumn(size_t start, size_t end) const -> ColumnView
        {
            return data_.subspan(start, end-start);
        }

    private:
        std::span<double> data_;
    };

    auto operator*(ColumnView lhs, ColumnView rhs) noexcept -> double;
    auto operator*(ColumnView lhs, double scalar) noexcept -> ColumnView;
}
