#include <iostream>
#include <vector>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "Percolation.h"

using namespace std;

int main()
{
	int n = 30;
	bool toDraw = true;

	srand(time(NULL));

	try {
		Percolation percolation(n);

		while(!percolation.percolates())
		{
			int i = n * (double)rand() / RAND_MAX;
			int j = n * (double)rand() / RAND_MAX;

			int row = i < n ? i : n - 1;
			int col = j < n ? j : n - 1;

			if (!percolation.isOpen(row, col))
			{
				percolation.open(row, col);
			}
		}

		if (toDraw)
		{
			cout << "'-': blocked site \n" << "' ': empty open site \n" << "'o': full open site\n" << endl;

			percolation.draw();
		}

		cout << "\np = " << (double)percolation.numberOfOpenSites() / (n * n) << endl;
	}
	catch (exception& ex) {
		cout << ex.what() << endl;
		exit(-1);
	}
	return 0;
}