#pragma once

#include <vector>
#include <string>

enum class Mode
{
    Linear,
    Quadratic,
    Cubic
};

class Interpolator
{
public:
    void LoadData(std::string filename);

    void Interpolate(Mode mode = Mode::Linear);

    std::vector<double> X();
    std::vector<double> Y();

private:
    // The number of intermediate points
    int _nInt = 20;

    // Raw data points
    std::vector<double> _xData;
    std::vector<double> _yData;

    // Interpolation result
    std::vector<double> _X;
    std::vector<double> _Y;
};