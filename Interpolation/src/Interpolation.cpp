#include <iostream>
#include <fstream>
#include <string>

#include "interpolationalgorithms.h"

using namespace std;

void Export(vector<double> x, vector<double> y, string filename)
{
    ofstream output;
    output.open(filename);
    for (int i = 0; i < x.size(); i++)
    {
        output << x[i] << " " << y[i] << endl;
    }
    output.close();
}

int main()
{
    Interpolator interpolator;

    interpolator.LoadData("../data/data.txt");

    interpolator.Interpolate(Mode::Linear);
    Export(interpolator.X(), interpolator.Y(), "linear.txt");

    interpolator.Interpolate(Mode::Quadratic);
    Export(interpolator.X(), interpolator.Y(), "quadratic.txt");

    interpolator.Interpolate(Mode::Cubic);
    Export(interpolator.X(), interpolator.Y(), "cubic.txt");
}
