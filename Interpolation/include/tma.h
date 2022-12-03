#pragma once

#include<vector>

// Tridiagonal Matrix Algorithm:
// M : matrix size
// a : lower diagonal
// b : main diagonal
// c : upper diagonal
// d : right-hand side

std::vector<double> TMA(
    size_t M,
    std::vector<double> a,
    std::vector<double> b,
    std::vector<double> c,
    std::vector<double> d);