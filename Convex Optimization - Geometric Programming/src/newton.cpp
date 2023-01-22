#include "newton.h"

#include "eigen-3.4.0/Eigen/Dense"

#include <iostream>
#include <chrono>

std::vector<double> Newton(std::vector<double> x, double (*func)(std::vector<double> x))
{
	bool log = false;

	int N = x.size();

	std::vector<double> step(N, 1e-6);

	std::vector<double> xPP, xPM, xMP, xMM;
	std::vector<double> xP, xM;

	Eigen::MatrixXd H(N, N);
	Eigen::VectorXd gradF(N);

	Eigen::VectorXd incr(N);

	int it = 0;
	
	auto start = std::chrono::high_resolution_clock::now();

	do
	{
		for (int row = 0; row < N; row++)
		{
			for (int col = row; col < N; col++)
			{
				double hRow = step[row];
				double hCol = step[col];

				xPP = x;
				xPP[row] += hRow;
				xPP[col] += hCol;

				xPM = x;
				xPM[row] += hRow;
				xPM[col] -= hCol;

				xMP = x;
				xMP[row] -= hRow;
				xMP[col] += hCol;

				xMM = x;
				xMM[row] -= hRow;
				xMM[col] -= hCol;

				H(row, col) = (func(xPP) - func(xPM) - func(xMP) + func(xMM)) / (hRow * hCol);
			}
		}

		for (int row = 0; row < N; row++)
		{
			for (int col = 0; col < row; col++)
			{
				H(row, col) = H(col, row);
			}
		}

		for (int row = 0; row < N; row++)
		{
			double hRow = step[row];

			xP = x;
			xP[row] += hRow;

			xM = x;
			xM[row] -= hRow;

			gradF(row) = (func(xP) - func(xM)) / (2 * hRow);
		}

		incr = -H.inverse() * gradF;

		for (int row = 0; row < N; row++)
			x[row] += incr(row);
		
		it++;

		if (log)
		{
			std::cout << it << " " 
				<< std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count() 
				<< " " << gradF.norm() << "\n";
		}
	}
	while (gradF.norm() > 1e-7);
	
	return x;
}