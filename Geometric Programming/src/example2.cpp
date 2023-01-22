#include <iostream>
#include <vector>

#include "newton.h"

using namespace std;

double Indicator(double val, double t)
{
	if (val >= 0)
		return 1e8;

	return -(1.0 / t) * log(-val);
}

double T;

double p = 1.0;

double ObjectiveFunc(vector<double> x)
{
	double val = 0;

	val += -(x[0] + x[1]);

	val += Indicator(
		log(
			(exp(x[0]) + exp(x[1])) * 2 / p
		), T);

	return val;
}

int main()
{
	vector<double> x = { log(0.1), log(0.1) };

	for (double t : {1, 4, 16, 64, 256, 1024, 4096})
	{
		T = t;
		x = Newton(x, &ObjectiveFunc);
	}
	
	cout << exp(x[0]) << " " << exp(x[1]) << "\n";
}