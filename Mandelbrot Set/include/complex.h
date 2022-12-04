#pragma once

class Complex
{
public:
    Complex();

    Complex(double x, double y);

    double Abs() const;

    Complex operator+(const Complex& b) const;

    Complex operator*(const Complex& b) const;

public:
    double _x;
    double _y;
};