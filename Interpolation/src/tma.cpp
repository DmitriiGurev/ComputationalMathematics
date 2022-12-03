#include "tma.h"

#include <stdexcept>

std::vector<double> TMA(
    size_t M,
    std::vector<double> a,
    std::vector<double> b,
    std::vector<double> c,
    std::vector<double> d)
{
    bool argsAreValid = M > 0 &&
        a.size() == M - 1 &&
        b.size() == M &&
        c.size() == M - 1 && 
        d.size() == M;

    if (!argsAreValid)
        throw std::invalid_argument("Invalid argument in TMA");

    std::vector<double> x(M);

    std::vector<double> p(M - 1);
    std::vector<double> q(M - 1);

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