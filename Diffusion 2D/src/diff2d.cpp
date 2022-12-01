#include "diff2d.h"

#include <vector>
#include <stdexcept>
#include <iostream>

Diffusion::Diffusion(double d, double tStep, UniformGrid grid, NeumannBC bc) :
	_d(d), _tStep(tStep), _grid(grid), _bc(bc)
{
	// Checking if the arguments are valid
	try
	{
		bool gridIsValid =
			_grid.lx > 0.0 &&
			_grid.ly > 0.0 &&
			_grid.nx > 0 &&
			_grid.ny > 0;

		if (!gridIsValid)
			throw std::invalid_argument("The grid parameters are required to be positive numbers.");
		if (d <= 0.0)
			throw std::invalid_argument("The diffusion coefficient is required to be a positive number.");
		if (tStep <= 0.0)
			throw std::invalid_argument("The time step is required to be a positive number.");
	}
	catch (std::invalid_argument& e)
	{
		std::cerr << e.what() << std::endl;
		exit(-1);
	}

	// Assembling the system using Triplets (line, row, value)
	double hx = grid.lx / (grid.nx - 1);
	double hy = grid.ly / (grid.ny - 1);

	double ax = d * tStep / (hx * hx);
	double ay = d * tStep / (hy * hy);

	typedef Eigen::Triplet<double> T;

	int nEq = grid.nx * grid.ny;

	_rhs = Eigen::VectorXd::Constant(nEq, 0.0);
	std::vector<T> coeffs;
	int id;

	// Edges: (0, 0)
	id = Id(0, 0);
	coeffs.push_back(T(id, Id(0, 0), 1 + 2 * ax + 2 * ay));

	coeffs.push_back(T(id, Id(1, 0), -ax));
	coeffs.push_back(T(id, Id(1, 0), -ax));
	_rhs(id) += -2 * ax * hx * bc.left;

	coeffs.push_back(T(id, Id(0, 1), -ay));
	coeffs.push_back(T(id, Id(0, 1), -ay));
	_rhs(id) += -2 * ay * hy * bc.bottom;

	// Edges: (0, ny - 1)
	id = Id(0, grid.ny - 1);
	coeffs.push_back(T(id, Id(0, grid.ny - 1), 1 + 2 * ax + 2 * ay));

	coeffs.push_back(T(id, Id(1, grid.ny - 1), -ax));
	coeffs.push_back(T(id, Id(1, grid.ny - 1), -ax));
	_rhs(id) += -2 * ax * hx * bc.left;

	coeffs.push_back(T(id, Id(0, grid.ny - 2), -ay));
	coeffs.push_back(T(id, Id(0, grid.ny - 2), -ay));
	_rhs(id) += 2 * ay * hy * bc.top;

	// Edges: (nx - 1, 0)
	id = Id(grid.nx - 1, 0);
	coeffs.push_back(T(id, Id(grid.nx - 1, 0), 1 + 2 * ax + 2 * ay));

	coeffs.push_back(T(id, Id(grid.nx - 2, 0), -ax));
	coeffs.push_back(T(id, Id(grid.nx - 2, 0), -ax));
	_rhs(id) += 2 * ax * hx * bc.right;

	coeffs.push_back(T(id, Id(grid.nx - 1, 1), -ay));
	coeffs.push_back(T(id, Id(grid.nx - 1, 1), -ay));
	_rhs(id) += -2 * ay * hy * bc.bottom;

	// Edges: (nx - 1, ny - 1)
	id = Id(grid.nx - 1, grid.ny - 1);
	coeffs.push_back(T(id, Id(grid.nx - 1, grid.ny - 1), 1 + 2 * ax + 2 * ay));

	coeffs.push_back(T(id, Id(grid.nx - 2, grid.ny - 1), -ax));
	coeffs.push_back(T(id, Id(grid.nx - 2, grid.ny - 1), -ax));
	_rhs(id) += 2 * ax * hx * bc.right;

	coeffs.push_back(T(id, Id(grid.nx - 1, grid.ny - 2), -ay));
	coeffs.push_back(T(id, Id(grid.nx - 1, grid.ny - 2), -ay));
	_rhs(id) += 2 * ay * hy * bc.top;

	// Walls: left
	for (int j = 1; j < grid.ny - 1; j++)
	{
		id = Id(0, j);
		coeffs.push_back(T(id, Id(0, j), 1 + 2 * ax + 2 * ay));

		coeffs.push_back(T(id, Id(1, j), -ax));
		coeffs.push_back(T(id, Id(1, j), -ax));
		_rhs(id) += -2 * ax * hx * bc.left;

		coeffs.push_back(T(id, Id(0, j + 1), -ay));
		coeffs.push_back(T(id, Id(0, j - 1), -ay));
	}

	// Walls: right
	for (int j = 1; j < grid.ny - 1; j++)
	{
		id = Id(grid.nx - 1, j);
		coeffs.push_back(T(id, Id(grid.nx - 1, j), 1 + 2 * ax + 2 * ay));

		coeffs.push_back(T(id, Id(grid.nx - 2, j), -ax));
		coeffs.push_back(T(id, Id(grid.nx - 2, j), -ax));
		_rhs(id) += 2 * ax * hx * bc.right;

		coeffs.push_back(T(id, Id(grid.nx - 1, j + 1), -ay));
		coeffs.push_back(T(id, Id(grid.nx - 1, j - 1), -ay));
	}

	// Walls: bottom
	for (int i = 1; i < grid.nx - 1; i++)
	{
		id = Id(i, 0);
		coeffs.push_back(T(id, Id(i, 0), 1 + 2 * ax + 2 * ay));

		coeffs.push_back(T(id, Id(i + 1, 0), -ax));
		coeffs.push_back(T(id, Id(i - 1, 0), -ax));

		coeffs.push_back(T(id, Id(i, 1), -ay));
		coeffs.push_back(T(id, Id(i, 1), -ay));
		_rhs(id) += -2 * ay * hy * bc.bottom;
	}

	// Walls: top
	for (int i = 1; i < grid.nx - 1; i++)
	{
		id = Id(i, grid.ny - 1);
		coeffs.push_back(T(id, Id(i, grid.ny - 1), 1 + 2 * ax + 2 * ay));

		coeffs.push_back(T(id, Id(i + 1, grid.ny - 1), -ax));
		coeffs.push_back(T(id, Id(i - 1, grid.ny - 1), -ax));

		coeffs.push_back(T(id, Id(i, grid.ny - 2), -ay));
		coeffs.push_back(T(id, Id(i, grid.ny - 2), -ay));
		_rhs(id) += 2 * ay * hy * bc.top;
	}

	// Non-boundary cells
	for (int i = 1; i < grid.nx - 1; i++)
	{
		for (int j = 1; j < grid.ny - 1; j++)
		{
			id = Id(i, j);
			coeffs.push_back(T(id, Id(i, j), 1 + 2 * ax + 2 * ay));
			coeffs.push_back(T(id, Id(i + 1, j), -ax));
			coeffs.push_back(T(id, Id(i - 1, j), -ax));
			coeffs.push_back(T(id, Id(i, j + 1), -ay));
			coeffs.push_back(T(id, Id(i, j - 1), -ay));
		}
	}

	_system = Eigen::SparseMatrix<double>(nEq, nEq);
	_system.setFromTriplets(coeffs.begin(), coeffs.end());
	_system.makeCompressed();

	// Setting up the solver
	_solver.analyzePattern(_system);
	_solver.factorize(_system);
}

std::vector<double> Diffusion::Solve(std::vector<double> u0, double tEnd)
{
	// Checking if the arguments are valid
	try
	{
		if (u0.size() != _grid.nx * _grid.ny)
			throw std::invalid_argument("The vector u0 must be of length nx * ny = " + std::to_string(_grid.nx * _grid.ny) + ".");
		if (tEnd <= 0.0)
			throw std::invalid_argument("The ending time is required to be a positive number.");
	}
	catch (std::invalid_argument& e)
	{
		std::cerr << e.what() << std::endl;
		exit(-1);
	}

	// Performing time iterations
	Eigen::VectorXd uPrev = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(u0.data(), u0.size());
	Eigen::VectorXd u(_grid.nx * _grid.ny);

	int nIt = tEnd / _tStep;

	for (int it = 0; it < nIt; it++)
	{
		u = _solver.solve(_rhs + uPrev);
		uPrev = u;
	}

	return std::vector<double>(u.data(), u.data() + u.rows() * u.cols());
}

int Diffusion::Id(int i, int j)
{
	return i + j * _grid.nx;
}