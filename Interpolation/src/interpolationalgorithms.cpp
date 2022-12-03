#include "interpolationalgorithms.h"
#include "tma.h"

#include <fstream>
#include <iostream>
#include <map>
#include <algorithm>

void Interpolator::LoadData(std::string filename)
{
    std::ifstream ifile(filename, std::ios::in);
    _xData.clear();
    _yData.clear();

    std::vector<double> xData;
    std::vector<double> yData;

    std::map<double, double> data;

    if (ifile.is_open())
    {
        double num = 0.0;
        while (ifile >> num)
        {
            xData.push_back(num);
            ifile >> num;
            yData.push_back(num);
        }

        transform(xData.begin(), xData.end(), yData.begin(), inserter(data, data.end()), [](double a, double b)
        {
            return std::make_pair(a, b);
        });

        for (auto element : data)
        {
            _xData.push_back(element.first);
            _yData.push_back(element.second);
        }
    }
}

void Interpolator::Interpolate(Mode mode)
{
    size_t N = _xData.size();

    if (N < 2)
    {
        std::cerr << "Too few data points for interpolation (N_min = 2);\n" <<
            "Interpolation failed.\n\n";
        return;
    }
    if (mode == Mode::Quadratic && N < 3)
    {
        mode = Mode::Linear;
        std::cerr << "Too few data points for quadratic interpolation (N_min = 3);\n" << 
            "The interpolation mode has been switched to Linear.\n\n";
    }
    if(mode == Mode::Cubic && N < 4)
    {
        mode = Mode::Linear;
        std::cerr << "Too few data points for cubic interpolation (N_min = 4);\n" <<
            "The interpolation mode has been switched to Linear.\n\n";
    }

    double h = _xData[1] - _xData[0]; 

    for (int i = 1; i + 1 < N; i++)
    {
        if (_xData[i + 1] - _xData[i] != h)
        {
            std::cerr << "The grid step must be constant;\n" <<
                "Interpolation failed.\n\n";
            return;
        }
    }

    _X.clear();
    _Y.clear();

    switch(mode)
    {
    case Mode::Linear:
        for (size_t i = 0; i < N - 1; i++)
        {
            for (size_t j = 0; j < _nInt; j++)
            {
                double x = _xData[0] + i * h + j * h / _nInt;
                double alpha = (double)j / _nInt;
                double y = (1 - alpha) * _yData[i] + alpha * _yData[i + 1];
                _X.push_back(x);
                _Y.push_back(y);
            }
        }
        break;

    case Mode::Quadratic:
        for (size_t i = 0; i + 2 < N; i += 2)
        {
            double a = _yData[i];
            double b = (_yData[i + 1] - _yData[i]) / h;
            double c = (_yData[i + 2] - 2 * _yData[i + 1] + _yData[i]) / (2 * h * h);

            for (size_t j = 0; j < _nInt; j++)
            {
                double x = _xData[0] + i * h + 2 * j * h / _nInt;
                double y = a + b * (x - _xData[i]) + c * (x - _xData[i]) * (x - _xData[i + 1]);
                _X.push_back(x);
                _Y.push_back(y);
            }
        }
        break;

    case Mode::Cubic:
        if (_xData.size() > 3)
        {
            std::vector<double> k_l(N - 3, h / 3);
            std::vector<double> k_c(N - 2, 4 * h / 3);
            std::vector<double> k_r(N - 3, h / 3);
            std::vector<double> RHS(N - 2);

            for (size_t i = 0; i + 2 < N; i++)
                RHS[i] = (_yData[i + 2] - 2 * _yData[i + 1] + _yData[i]) / h;

            std::vector<double> c_ = TMA(N - 2, k_l, k_c, k_r, RHS);
            std::vector<double> c(N);

            c[0] = 0;
            for (size_t i = 1; i + 1 < N; i++)
                c[i] = c_[i - 1];
            c[N - 1] = 0;

            std::vector<double> d(N);
            for (size_t i = 0; i + 1 < N; i++)
                d[i] = (c[i + 1] - c[i]) / (3 * h);

            std::vector<double> b(N);
            for (size_t i = 0; i + 1 < N; i++)
                b[i] = (_yData[i + 1] - _yData[i]) / h - (h / 3) * (c[i + 1] + 2 * c[i]);

            for (size_t i = 0; i + 1 < N; i++)
            {
                for (size_t j = 0; j < _nInt; j++)
                {
                    double x = _xData[0] + i * h + j * h / _nInt;
                    double y = _yData[i] + b[i] * (x - _xData[i]) + c[i] * (x - _xData[i]) * (x - _xData[i]) +
                        d[i] * (x - _xData[i]) * (x - _xData[i]) * (x - _xData[i]);
                    _X.push_back(x);
                    _Y.push_back(y);
                }
            }
        }
        break;
    }
}

std::vector<double> Interpolator::X()
{
    return _X;
}

std::vector<double> Interpolator::Y()
{
    return _Y;
}