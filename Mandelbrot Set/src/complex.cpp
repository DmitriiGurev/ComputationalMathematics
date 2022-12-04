#include "complex.h"

#include <cmath>

Complex::Complex() :
    _x(0), _y(0)
{}

Complex::Complex(double x, double y) :
    _x(x), _y(y)
{}

double Complex::Abs() const
{
    return sqrt(_x * _x + _y * _y);
}

Complex Complex::operator+(const Complex& b) const
{
    return Complex(_x + b._x, _y + b._y);
}

Complex Complex::operator*(const Complex& b) const
{
    return Complex(_x * b._x - _y * b._y, _y * b._x + _x * b._y);
}