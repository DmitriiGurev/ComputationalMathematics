#include <solver.h>

#include <iostream>
#include <fstream>

using namespace std;

int main()
{
    Solver solver;
    solver.xN = 1e3;
    solver.xStep = 2 * 1.0e-3;

    double rhoBase = 1.0e-3;

    InitConditions initCond;
    initCond.rho = vector<double>(solver.xN, rhoBase);
    initCond.u = vector<double>(solver.xN, 0);
    for (int i = 0; i < solver.xN; i++)
    {
        double x = solver.xStep * i;
        if (x <= 1)
            initCond.rho[i] = 1;

        if (x <= 0.4)
            initCond.rho[i] = 0.5;
    }

    BndConditions bndCond;
    bndCond.uLeft = 0;
    bndCond.uRight = 0;
    bndCond.rhoRight = rhoBase;

    solver.tStep = 2.5 * 1.0e-4;
    solver.tN = 8.5 / solver.tStep;

    solver.writingStep = 50;

    Solution solution = solver.Solve<WenoLF>(initCond, bndCond);
}