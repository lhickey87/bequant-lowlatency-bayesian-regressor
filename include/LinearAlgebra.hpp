#pragma once
#include <vector>
#include <cassert>

namespace LA {
    void partialCopy(std::vector<double>& data, std::vector<double>& result, size_t length, size_t index);
    auto dotProduct(std::vector<double>& lhs, std::vector<double>& rhs);
    auto partialDotProduct(std::vector<double>& lhs, std::vector<double>& rhs,size_t length, size_t index);
}
