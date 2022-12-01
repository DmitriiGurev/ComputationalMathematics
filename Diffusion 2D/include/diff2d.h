#pragma once

#include <eigen-3.4.0/Eigen/Sparse>

struct UniformGrid
{
	double lx = 0.0;
	double ly = 0.0;
	int	   nx = 0;
	int	   ny = 0;
};

struct NeumannBC
{
	double left = 0.0;
	double bottom = 0.0;
	double right = 0.0;
	double top = 0.0;
};

class Diffusion
{
public:
	Diffusion(double d,	double tStep, UniformGrid grid, NeumannBC bc);

	std::vector<double> Solve(std::vector<double> u0, double tEnd);

private:
	// Index of the (i,j) grid point
	int Id(int i, int j);

private:
	double		_d;
	double		_tStep;
	UniformGrid _grid;
	NeumannBC	_bc;

	Eigen::VectorXd				_rhs;
	Eigen::SparseMatrix<double> _system;

	Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> _solver;
};