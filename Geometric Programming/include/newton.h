#pragma once

#include <vector>

std::vector<double> Newton(std::vector<double> x, double (*func)(std::vector<double> x));