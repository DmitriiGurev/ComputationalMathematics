#pragma once

#include <vector>
#include <string>

enum Scheme
{
    LaxFriedrichs,
    LaxWendroff,
    LargeParticle,
    WenoLF
};

struct Solution
{
    std::vector<double> rho;
    std::vector<double> u;
    std::vector<double> p;
};

struct InitConditions
{
    std::vector<double> rho;
    std::vector<double> u;
};

struct BndConditions
{
    double uLeft;
    double uRight;
    double rhoRight;
};

class Solver
{
public:
    Solver() {}

    Solver(int xN, double xStep, int tN, double tStep) :
        xN(xN), xStep(xStep), tN(tN), tStep(tStep) {}

    template <Scheme scheme>
    Solution Solve(const InitConditions& intiCond, const BndConditions& bndCond);

    void WriteResults(std::string postfix, const Solution& sol);

    void WriteToVTK(std::string fileName, const Solution& sol);

    void Smooth(Solution& sol, double a);

public:
    int xN;
    int tN;

    double xStep;
    double tStep;

    double gamma = 5.0 / 3;

    int writingStep = 5000;

    int smoothingStep = 1;
    double smoothingAlpha = 0.01;
};