#pragma once
#include <span>
#include <numeric>
#include <cassert>

namespace LA {
    template <typename T>
    struct ColumnView {

        using span_t = std::span<T>;

        ColumnView(span_t col) : data_(col) {}
        explicit ColumnView(T* data, size_t length) : data_(data, length) {}

        auto inline operator[](size_t index) -> T&
        {
            return data_[index];
        }

        auto inline operator[](size_t index) const -> const T&
        {
            return data_[index];
        }

        auto operator-=(ColumnView rhs) -> void
        {
            assert(length() == rhs.length());
            for (size_t i = 0; i < length(); ++i)
            {
                data_[i] -= rhs[i];
            }
        }

        auto operator*=(double scalar) -> void
        {
            for (size_t i = 0; i < length();++i)
            {
                data_[i] *= scalar;
            }
        }

        auto inline begin() {return data_.begin();}
        auto inline end() {return data_.end();}
        auto inline begin() const {return data_.begin();}
        auto inline end() const {return data_.end();}
        auto inline length() const -> size_t{ return data_.size();}

        auto inline data() const -> T* { return data_.data(); }

        auto inline addScaled(ColumnView x, double a) noexcept
        {
            assert(length() == x.length());
            for (size_t i = 0; i < length(); ++i) {
                data_[i] += a * x[i];
            }
        }

        auto inline SubColumn(size_t start, size_t end) const -> ColumnView
        {
            return data_.subspan(start, end-start);
        }

    private:
        span_t data_;
    };

    using Column = ColumnView<double>;
    using ConstColumn = ColumnView<const double>;

    template <typename T>
    auto operator*(ColumnView<T> lhs, ColumnView<T> rhs) noexcept -> double
    {
        double sum = 0.0;
        for (size_t i = 0; i < lhs.length();++i)
        {
            sum += lhs[i]*rhs[i];
        }
        return sum;
    }

    inline auto operator*(Column lhs, double scalar) noexcept -> Column
    {
        for (size_t i = 0; i < lhs.length(); ++i)
        {
            lhs[i] *= scalar;
        }
        return lhs;
    }
}
