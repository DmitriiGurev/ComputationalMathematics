#include <stdexcept>
#include <iostream>

#include "Percolation.h"

int Percolation::Index(int row, int col)
{
	return row * _n + col;
}

Percolation::Percolation(int n) : _UF(n * n + 2)
{
	if (n < 0)
		throw std::runtime_error("n must be more than 0");
	_n = n;

	// Sites: [0] to [n^2 - 1]
	// Virtual sites: [n^2] - top row; [n^2 + 1] - bottom row
	_open.resize(_n * _n, false);
}

void Percolation::open(int row, int col)
{
	if (row >= _n || _n < 0 || col >= _n || col < 0)
		throw std::runtime_error("row, col must be in [0, N-1]");
	if (isOpen(row, col)) return;
	int index = Index(row, col);
	_open[index] = true;
	if (row == 0 && col == 0)
	{
		if (isOpen(0, 1)) _UF.Union(index, Index(0, 1));
		if (isOpen(1, 0)) _UF.Union(index, Index(1, 0));
		_UF.Union(index, _n * _n);
	}
	else if (row == 0 && col == _n - 1)
	{
		if (isOpen(0, _n - 2)) _UF.Union(index, Index(0, _n - 2));
		if (isOpen(1, _n - 1)) _UF.Union(index, Index(1, _n - 1));
		_UF.Union(index, _n * _n);
	}
	else if (row == _n - 1 && col == 0)
	{
		if (isOpen(_n - 2, 0)) _UF.Union(index, Index(_n - 2, 0));
		if (isOpen(_n - 1, 1)) _UF.Union(index, Index(_n - 1, 1));
		_UF.Union(index, _n * _n + 1);
	}
	else if (row == _n - 1 && col == _n - 1)
	{
		if (isOpen(_n - 1, _n - 2)) _UF.Union(index, Index(_n - 1, _n - 2));
		if (isOpen(_n - 2, _n - 1)) _UF.Union(index, Index(_n - 2, _n - 1));
		_UF.Union(index, _n * _n + 1);
	}
	else if (row == 0 && col > 0 && col < _n - 1)
	{
		if (isOpen(0, col - 1)) _UF.Union(index, Index(0, col - 1));
		if (isOpen(0, col + 1)) _UF.Union(index, Index(0, col + 1));
		if (isOpen(1, col)) _UF.Union(index, Index(1, col));
		_UF.Union(index, _n * _n);
	}
	else if (row == _n - 1 && col > 0 && col < _n - 1)
	{
		if (isOpen(_n - 1, col - 1)) _UF.Union(index, Index(_n - 1, col - 1));
		if (isOpen(_n - 1, col + 1)) _UF.Union(index, Index(_n - 1, col + 1));
		if (isOpen(_n - 2, col)) _UF.Union(index, Index(_n - 2, col));
		_UF.Union(index, _n * _n + 1);
	}
	else if (row > 0 && row < _n - 1 && col == 0)
	{
		if (isOpen(row - 1, 0)) _UF.Union(index, Index(row - 1, 0));
		if (isOpen(row + 1, 0)) _UF.Union(index, Index(row + 1, 0));
		if (isOpen(row, 1)) _UF.Union(index, Index(row, 1));
	}
	else if (row > 0 && row < _n - 1 && col == _n - 1)
	{
		if (isOpen(row - 1, _n - 1)) _UF.Union(index, Index(row - 1, _n - 1));
		if (isOpen(row + 1, _n - 1)) _UF.Union(index, Index(row + 1, _n - 1));
		if (isOpen(row, _n - 2)) _UF.Union(index, Index(row, _n - 2));
	}
	else
	{
		if (isOpen(row, col + 1)) _UF.Union(index, Index(row, col + 1));
		if (isOpen(row, col - 1)) _UF.Union(index, Index(row, col - 1));
		if (isOpen(row + 1, col)) _UF.Union(index, Index(row + 1, col));
		if (isOpen(row - 1, col)) _UF.Union(index, Index(row - 1, col));
	}
}

bool Percolation::isOpen(int row, int col)
{
	return _open[Index(row, col)];
}

bool Percolation::isFull(int row, int col)
{
	return _UF.Connected(Index(row, col), _n * _n);
}

int Percolation::numberOfOpenSites()
{
	return count(begin(_open), end(_open), true);
}

bool Percolation::percolates()
{
	return _UF.Connected(_n * _n, _n * _n + 1);
}

void Percolation::draw()
{
	for (int i = 0; i < _n; i++)
	{
		for (int j = 0; j < _n; j++)
		{
			if (isOpen(i, j) && isFull(i, j))
			{
				std::cout << "o" << " ";
			}
			else if (isOpen(i, j))
			{
				std::cout << " " << " ";
			}
			else
			{
				std::cout << "-" << " ";
			}
		}
		std::cout << "\n";
	}
}