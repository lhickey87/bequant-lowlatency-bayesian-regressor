#pragma once
#include "Matrix.hpp"
#include <iostream>

class OLS {

    auto fit() noexcept;
    auto solve() noexcept;

public:

private:
    std::vector<double> betas;
    size_t size;
};
