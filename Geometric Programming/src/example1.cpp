#include <iostream>
#include <vector>

#include "newton.h"

#include "eigen-3.4.0/Eigen/Dense"

using namespace std;

Eigen::Matrix4d M;

double ObjectiveFunc(vector<double> x)
{
	double val = 0;

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			val += exp(2 * x[i] - 2 * x[j] + log(M(i, j) * M(i, j)));
		}
	}
	return log(val);
}

int main()
{
	vector<double> x(4, -1);

	M << 1, 2, 3, 4,
		5, 6, 7, 8,
		9, 10, 11, 12,
		13, 14, 15, 16;

	cout << "M =\n" << M << "\n\n";

	x = Newton(x, &ObjectiveFunc);

	Eigen::MatrixXd D = Eigen::MatrixXd::Zero(4, 4);
	for (int i = 0; i < 4; i++)
		D(i, i) = exp(x[i]);

	cout << "D =\n" << D << "\n\n";

	cout << "||M|| = " << M.norm() << "\n\n";

	cout << "||D M D^-1|| = " << (D * M * D.inverse()).norm();
}