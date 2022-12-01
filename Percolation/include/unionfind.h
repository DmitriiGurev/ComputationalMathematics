#pragma once

#include <vector>

class UnionFind
{
// Union-Find: Weighted Quick-Union with Path Compression (WQUPC)
public:
	UnionFind(int N);

	bool Connected(int p, int q);

	void Union(int p, int q);

	void Print() const;

private:
	std::vector<int> id;
	std::vector<int> sz;
	int Root(int i);
};