#include <iostream>
#include <vector>
#include <fstream>

#include "newton.h"

using namespace std;

double Indicator(double val, double t)
{
	if (val >= 0)
		return 1e8;

	return -(1.0 / t) * log(-val);
}

double wMin, wMax;
double hMin, hMax;

double F;
double E;
double sigmaMax;
double yMax;

double T;

double ObjectiveFunc(vector<double> x)
{
	double val = 0;

	val += log(exp(x[0] + x[3]) + exp(x[1] + x[4]) + exp(x[2] + x[5]));

	for (int i : {0, 1, 2})
	{
		val += Indicator(x[i] - log(wMax), T);
		val += Indicator(log(wMin) - x[i], T);
		val += Indicator(x[i + 3] - log(hMax), T);
		val += Indicator(log(hMin) - x[i + 3], T);

		val += Indicator(-x[i] - 2 * x[i + 3] + log(6 * (i + 1) * F / sigmaMax), T);
	}
	
	val += Indicator(
		log(
			exp(-x[0] - 3 * x[3] + log(4 * F / (E * yMax)))
			+ exp(-x[1] - 3 * x[4] + log(28 * F / (E * yMax)))
			+ exp(-x[2] - 3 * x[5] + log(76 * F / (E * yMax)))
		), T);
	
	return val;
}

int main()
{
	wMin = 0.1;			// [m]
	wMax = 0.5;			// [m]
	hMin = 0.1;			// [m]
	hMax = 0.5;			// [m]

	F = 3e5;			// [N]
	E = 200e9;			// [Pa]
	sigmaMax = 100e6;	// [Pa] 
	yMax = 0.01;		// [m]

	vector<double> x(6, 0.49);

	for (int i = 0; i < 6; i++)
		x[i] = log(x[i]);

	for (double t : { 1, 4, 16, 64, 256, 1024, 4096 })
	{
		T = t;
		x = Newton(x, &ObjectiveFunc);
	}

	double w1 = exp(x[0]);
	double w2 = exp(x[1]);
	double w3 = exp(x[2]);
	double h1 = exp(x[3]);
	double h2 = exp(x[4]);
	double h3 = exp(x[5]);

	cout << "w1 = " << w1 << " m\n";
	cout << "h1 = " << h1 << " m\n";
	cout << "w2 = " << w2 << " m\n";
	cout << "h2 = " << h2 << " m\n";
	cout << "w3 = " << w3 << " m\n";
	cout << "h3 = " << h3<< " m\n\n";

	cout << "S = " << w1 * h1 + w2 * h2 + w3 * h3 << " m^2\n\n";

	cout << "s1 / sMax = " << 6 * 1 * F / (w1 * h1 * h1 * sigmaMax) << "\n";
	cout << "s2 / sMax = " << 6 * 2 * F / (w2 * h2 * h2 * sigmaMax) << "\n";
	cout << "s3 / sMax = " << 6 * 3 * F / (w3 * h3 * h3 * sigmaMax) << "\n\n";

	cout << "y1 / yMax = " << 4 * F / (E * w1 * h1 * h1 * h1 * yMax) + 
					   28 * F / (E * w2 * h2 * h2 * h2 * yMax) +
					   76 * F / (E * w3 * h3 * h3 * h3 * yMax);
}