#include <iostream>

#include "complex.h"

using namespace std;

Complex MandelbrotF(const Complex& z, const Complex& c)
{
    return z * z + c;
}

bool Diverges(const Complex& c, int nIt, double eps)
{
    Complex z(0, 0);

    for (int i = 0; i < nIt; i++)
    {
        z = MandelbrotF(z, c);
        if (z.Abs() > eps)
            return true;
    }
    return false;
}

int main()
{    
    int nIt = 3000;
    double eps = 100;

    double xStep = 4*0.00005, yStep = 4*0.0001;
    double xMin = 0.3675, xMax = 0.420;
    double yMin = 0.188, yMax = 0.239;

    int nX = (xMax - xMin) / xStep;
    int nY = (yMax - yMin) / yStep;

    for (int i = nY; i > -1; i--)
    {
        for (int j = 0; j < nX; j++)
        {
            Complex point(xMin + j * (xMax - xMin) / nX,
                yMin + i * (yMax - yMin) / nY);
            if (!Diverges(point, nIt, eps))
            {
                cout << "X";
            }
            else
            {
                cout << " ";
            }
        }
        cout << "\n";
    }
}