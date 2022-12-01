#pragma once

#include "UnionFind.h"

#include <vector>

class Percolation
{
public:
	Percolation(int n);

	// Opens the site (row, col) if it is not open already
	void open(int row, int col);

	// Is the site (row, col) open?
	bool isOpen(int row, int col);

	// Is the site (row, col) full?
	bool isFull(int row, int col);

	// Returns the number of open sites
	int numberOfOpenSites();

	// Does the system percolate?
	bool percolates();

	void draw();

private:
	int Index(int row, int col);

private:
	UnionFind _UF;
	std::vector<bool> _open;
	int _n;
};