#include "UnionFind.h"

#include <iostream>

UnionFind::UnionFind(int N)
{
	id.resize(N);
	sz.resize(N);
	for (int i = 0; i < N; i++)
	{
		id[i] = i;
		sz[i] = 1;
	}
}

bool UnionFind::Connected(int p, int q)
{
	return Root(p) == Root(q);
}

void UnionFind::Union(int p, int q)
{
	int i = Root(p);
	int j = Root(q);

	if (i == j) return;

	if (sz[i] < sz[j])
	{
		id[i] = j; sz[j] += sz[i];
	}
	else
	{
		id[j] = i; sz[i] += sz[j];
	}
}

void UnionFind::Print() const
{
	for (int i : id)
		std::cout << i << " ";
	std::cout << std::endl;
}

int UnionFind::Root(int i)
{
	while (id[i] != i)
	{
		id[i] = id[id[i]];
		i = id[i];
	}
	return i;
}