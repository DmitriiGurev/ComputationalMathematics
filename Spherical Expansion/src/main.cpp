#include <solver.h>

#include <iostream>
#include <fstream>

using namespace std;

int main()
{
    Solver solver;
    solver.xN = 1e4;
    solver.xStep = 1.0e-3;

    InitConditions initCond;
    initCond.rho = vector<double>(solver.xN, 1.0e-3);
    initCond.u = vector<double>(solver.xN, 0);
    for (int i = 0; i < solver.xN; i++)
    {
        double x = solver.xStep * i;
        if (x <= 1)
            initCond.rho[i] = 1;
    }

    BndConditions bndCond;
    bndCond.uLeft = 0;
    bndCond.uRight = 0;
    bndCond.rhoRight = 0;

    solver.smoothingStep = 1;
    solver.smoothingAlpha = 0.1;

    solver.tStep = 1.0e-4;
    solver.tN = 8.5 / solver.tStep;

    solver.writingStep = 1000;

    Solution solution = solver.Solve<LaxWendroff>(initCond, bndCond);
}