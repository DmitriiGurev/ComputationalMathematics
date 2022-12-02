#pragma once

#include <fstream>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <string>
#include <map>
#include <algorithm>

class Interpolator
{
public:
    Interpolator();

    void LoadData(std::string filename);

    virtual std::vector<double> GetOutput_x() = 0;
    virtual std::vector<double> GetOutput_y() = 0;

protected:
    virtual void update() = 0;

protected:
    std::vector<double> data_x;
    std::vector<double> data_y;
};

class Linear : public Interpolator
{
public:
    std::vector<double> GetOutput_x();
    std::vector<double> GetOutput_y();

private:
    void update();
};


class Quadratic : public Interpolator
{
public:
    std::vector<double> GetOutput_x();
    std::vector<double> GetOutput_y();

private:
    void update();

private:
    double n_submesh = 30;
    std::vector<double> output_x;
    std::vector<double> output_y;
};

class Cubic : public Interpolator
{
public:
    std::vector<double> GetOutput_x();
    std::vector<double> GetOutput_y();

private:
    void update();

private:
    double n_submesh = 20;
    std::vector<double> output_x;
    std::vector<double> output_y;
};