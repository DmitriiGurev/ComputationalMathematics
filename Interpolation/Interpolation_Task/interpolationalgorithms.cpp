#include "InterpolationAlgorithms.h"

using namespace std;

Interpolator::Interpolator() {}

void Interpolator::LoadData(string filename)
{
    ifstream ifile(filename, ios::in);
    data_x.clear();
    data_y.clear();

    vector<double> _data_x;
    vector<double> _data_y;

    map<double, double> data;

    if (ifile.is_open()) {
        double num = 0.0;
        while (ifile >> num) {
            _data_x.push_back(num);
            ifile >> num;
            _data_y.push_back(num);
        }

        transform(_data_x.begin(), _data_x.end(), _data_y.begin(), inserter(data, data.end()), [](double a, double b)
        {
            return std::make_pair(a, b);
        });

        for (auto element : data)
        {
            data_x.push_back(element.first);
            data_y.push_back(element.second);
        }
    }

    update();
}


// Linear
void Linear::update() {}

vector<double> Linear::GetOutput_x()
{
    return data_x;
}
vector<double> Linear::GetOutput_y()
{
    return data_y;
}

// Quadratic
void Quadratic::update()
{
    output_x.clear();
    output_y.clear();

    if (data_x.size() >= 3)
    {
        double h = data_x[1] - data_x[0]; // Assumed to be constant
        for (size_t i = 0; i + 2 < data_x.size(); i += 2)
        {
            double a = data_y[i];
            double b = (data_y[i + 1] - data_y[i]) / h;
            double c = (data_y[i + 2] - 2 * data_y[i + 1] + data_y[i]) / (2 * h * h);

            for (size_t j = 0; j < n_submesh; j++)
            {
                double x = data_x[0] + i * h + 2 * j * h / n_submesh;
                double y = a + b * (x - data_x[i]) + c * (x - data_x[i]) * (x - data_x[i + 1]);
                output_x.push_back(x);
                output_y.push_back(y);
            }
        }
    }
}

vector<double> Quadratic::GetOutput_x()
{
    return output_x;
}

vector<double> Quadratic::GetOutput_y()
{
    return output_y;
}

// Tridiagonal matrix algorithm
vector<double> TMA(size_t M, vector<double> a, vector<double> b, vector<double> c, vector<double> d)
{
    vector<double> x(M);

    vector<double> p(M - 1);
    vector<double> q(M - 1);

    p[0] = -c[0] / b[0];
    q[0] = d[0] / b[0];

    for (size_t i = 0; i < M - 2; i++)
    {
        p[i + 1] = -c[i + 1] / (a[i + 1] * p[i] + b[i + 1]);
        q[i + 1] = (d[i + 1] - a[i + 1] * q[i]) / (a[i + 1] * p[i] + b[i + 1]);
    }

    x[M - 1] = (d[M - 1] - a[M - 2] * q[M - 2]) / (p[M - 2] * a[M - 2] + b[M - 1]);

    for (size_t i = 0; i < M - 1; i++)
    {
        x[M - 2 - i] = p[M - 2 - i] * x[M - 1 - i] + q[M - 2 - i];
    }

    return x;
}

// Cubic
void Cubic::update() {
    output_x.clear();
    output_y.clear();

    if (data_x.size() > 3)
    {
        double h = data_x[1] - data_x[0]; // Assumed to be constant

        size_t N = data_x.size();
        vector<double> k_l(N - 3, h / 3);
        vector<double> k_c(N - 2, 4 * h / 3);
        vector<double> k_r(N - 3, h / 3);
        vector<double> RHS(N - 2);

        for (size_t i = 0; i < N - 2; i++)
            RHS[i] = (data_y[i + 2] - 2 * data_y[i + 1] + data_y[i]) / h;

        vector<double> c_ = TMA(N - 2, k_l, k_c, k_r, RHS); // l, r ?
        vector<double> c(N);

        c[0] = 0;
        for (size_t i = 1; i < N - 1; i++)
            c[i] = c_[i - 1];
        c[N - 1] = 0;

        vector<double> d(N);
        for (size_t i = 0; i < N - 1; i++)
            d[i] = (c[i + 1] - c[i]) / (3 * h);

        vector<double> b(N);
        for (size_t i = 0; i < N - 1; i++)
            b[i] = (data_y[i + 1] - data_y[i]) / h - (h / 3) * (c[i + 1] + 2 * c[i]);

        for (size_t i = 0; i < N - 1; i++)
        {
            for (size_t j = 0; j < n_submesh; j++)
            {
                double x = data_x[0] + i * h + j * h / n_submesh;
                double y = data_y[i] + b[i] * (x - data_x[i]) + c[i] * (x - data_x[i]) * (x - data_x[i]) +
                    d[i] * (x - data_x[i]) * (x - data_x[i]) * (x - data_x[i]);
                output_x.push_back(x);
                output_y.push_back(y);
            }
        }
    }
}

vector<double> Cubic::GetOutput_x()
{
    return output_x;
}

vector<double> Cubic::GetOutput_y()
{
    return output_y;
}