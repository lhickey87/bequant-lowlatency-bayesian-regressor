#include "LinearAlgebra.hpp"

namespace LA {

    void partialCopy(std::vector<double> &data, std::vector<double> &result, size_t length, size_t index){
        assert(data.size() == result.size());
        for (size_t i = 0; i < length; ++i) result[i] = data[i+index];
    }

    auto dotProduct(std::vector<double> &lhs, std::vector<double> &rhs)
    {
        assert(lhs.size() == rhs.size());
        double result = 0.0;
        for (int i = 0; i < lhs.size(); ++i)
        {
            result += lhs[i]*rhs[i];
        }
        return result;
    }

    auto partialDotProduct(std::vector<double> &lhs, std::vector<double> &rhs, size_t length,size_t index)
    {
        assert(lhs.size() == rhs.size());
        double result = 0.0;
        for (int i = index; i < length; ++i){
            result += lhs[i]*rhs[i];
        }

    }


}
