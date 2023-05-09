#include <solver.h>

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>

using namespace std;

void Solver::WriteResults(std::string postfix, const Solution& sol)
{
    ofstream out;
    out.open("solution/rho/rho_" + postfix + ".txt");
    if (out.is_open()) {
        for (int i = 0; i < xN; i++)
            out << i * xStep << " " << sol.rho[i] << "\n";
    }
    out.close();

    out.open("solution/u/u_" + postfix + ".txt");
    if (out.is_open()) {
        for (int i = 0; i < xN; i++)
            out << i * xStep << " " << sol.u[i] << "\n";
    }
    out.close();
    
    out.open("solution/p/p_" + postfix + ".txt");
    if (out.is_open()) {
        for (int i = 0; i < xN; i++)
            out << i * xStep << " " << sol.p[i] << "\n";
    }
    out.close();

    out.open("rho.txt");
    for (int i = 0; i < xN; i++)
        out << i * xStep << " " << sol.rho[i] << "\n";
    out.close();
}


void Solver::Smooth(Solution& sol, double a)
{
    Solution solSmooth = sol;
    for (int i = 1; i < xN - 1; i++)
    {
        solSmooth.u[i] = (1 - 2 * a) * sol.u[i] + a * sol.u[i + 1] + a * sol.u[i - 1];
        solSmooth.rho[i] = (1 - 2 * a) * sol.rho[i] + a * sol.rho[i + 1] + a * sol.rho[i - 1];
    }
    sol = solSmooth;
}

template<>
Solution Solver::Solve<LaxFriedrichs>(const InitConditions& intiCond, const BndConditions& bndCond)
{
    Solution sol;
    sol.u = intiCond.u;
    sol.rho = intiCond.rho;
    sol.p = vector<double>(xN);
    for (int i = 0; i < xN; i++)
        sol.p[i] = pow(sol.rho[i], gamma);


    for (int it = 0; it < tN; it++)
    {
        if (it % writingStep == 0)
        {
            cout << "t = " << it * tStep << "\n";
            double maxU = *max_element(sol.u.begin(), sol.u.end());
            cout << "CFL = " << tStep * maxU / xStep << "\n\n";

            WriteResults(to_string(it * tStep), sol);
        }

        Solution newSol;
        newSol.rho = vector<double>(xN);
        newSol.u = vector<double>(xN);
        newSol.p = vector<double>(xN);

        if (it % smoothingStep == 0)
            Smooth(sol, smoothingAlpha);

        for (int i = 1; i < xN - 1; i++)
        {
            newSol.u[i] = (sol.u[i + 1] + sol.u[i - 1]) / 2 - 
                tStep / (2 * xStep) * ((sol.p[i + 1] - sol.p[i - 1]) / 
                ((sol.rho[i + 1] + sol.rho[i - 1]) / 2) +
                ((sol.u[i + 1] + sol.u[i - 1]) / 2) * (sol.u[i + 1] - sol.u[i - 1]));

            if (newSol.u[i] != newSol.u[i])
                throw runtime_error("nan appeared");
        }
        newSol.u[0] = bndCond.uLeft;
        newSol.u[xN - 1] =  bndCond.uRight;

        double x;
        for (int i = 1; i < xN - 1; i++)
        {
            x = i * xStep;
            newSol.rho[i] = (sol.rho[i + 1] + sol.rho[i - 1]) / 2 - tStep / (2 * xStep) *
                (pow(1 + xStep / x, 2) * sol.rho[i + 1] * sol.u[i + 1] -
                pow(1 - xStep / x, 2) * sol.rho[i - 1] * sol.u[i - 1]);
        }      
        newSol.rho[0] = newSol.rho[1];
        newSol.rho[xN - 1] = bndCond.rhoRight;

        for (int i = 1; i < xN - 1; i++)
        {
            newSol.p[i] = pow((sol.rho[i + 1] + sol.rho[i - 1]) / 2, gamma);
        }
        newSol.p[0] = pow(sol.rho[0], gamma);
        newSol.p[xN - 1] = pow(sol.rho[xN - 1], gamma);

        sol = newSol;
    }
    return sol;
}

template<>
Solution Solver::Solve<LaxWendroff>(const InitConditions& intiCond, const BndConditions& bndCond)
{
    Solution sol;
    sol.u = intiCond.u;
    sol.rho = intiCond.rho;
    sol.p = vector<double>(xN);
    for (int i = 0; i < xN; i++)
        sol.p[i] = pow(sol.rho[i], gamma);

    double s = tStep / xStep; 

    for (int it = 0; it < tN; it++)
    {
        if (it % writingStep == 0)
        {
            cout << "t = " << it * tStep << "\n";
            double maxU = *max_element(sol.u.begin(), sol.u.end());
            cout << "CFL = " << tStep * maxU / xStep << "\n\n";

            WriteResults(to_string(it * tStep), sol);
        }

        Solution newSol;
        newSol.rho = vector<double>(xN);
        newSol.u = vector<double>(xN);
        newSol.p = vector<double>(xN);

        if (it % smoothingStep == 0)
            Smooth(sol, smoothingAlpha);

        for (int i = 1; i < xN - 1; i++)
        {
            double r = xStep * i;

            // TODO: Precompute these
            double uHalfUp = 0.5 * (sol.u[i + 1] + sol.u[i]) -
                0.25 * s * (sol.u[i + 1] + sol.u[i]) * (sol.u[i + 1] - sol.u[i]) - 
                0.5 * s * (sol.p[i + 1] - sol.p[i]) / (0.5 * (sol.rho[i + 1] + sol.rho[i]));

            double uHalfDown = 0.5 * (sol.u[i] + sol.u[i - 1]) -
                0.25 * s * (sol.u[i] + sol.u[i - 1]) * (sol.u[i] - sol.u[i - 1]) -
                0.5 * s * (sol.p[i] - sol.p[i - 1]) / (0.5 * (sol.rho[i] + sol.rho[i - 1]));

            double rhoHalfUp = 0.5 * (sol.rho[i + 1] + sol.rho[i]) -
                0.5 * s * (pow(r + xStep, 2) * sol.rho[i + 1] * sol.u[i + 1] - 
                pow(r, 2) * sol.rho[i] * sol.u[i]) / pow(r + 0.5 * xStep, 2);

            double rhoHalfDown = 0.5 * (sol.rho[i] + sol.rho[i - 1]) -
                0.5 * s * (pow(r, 2) * sol.rho[i] * sol.u[i] - 
                pow(r - xStep, 2) * sol.rho[i - 1] * sol.u[i - 1]) / pow(r - 0.5 * xStep, 2);

            newSol.u[i] = sol.u[i] - s * sol.u[i] * (uHalfUp - uHalfDown) - 
                s * gamma * pow(sol.rho[i], gamma - 2) * (rhoHalfUp - rhoHalfDown);

            newSol.rho[i] = sol.rho[i] - s * (pow(r + 0.5 * xStep, 2) * rhoHalfUp * uHalfUp - 
                pow(r - 0.5 * xStep, 2) * rhoHalfDown * uHalfDown) / pow(r, 2);

            if (newSol.u[i] != newSol.u[i])
                throw runtime_error("nan appeared");
        }

        newSol.u[0] = bndCond.uLeft;
        newSol.u[xN - 1] = bndCond.uRight;

        newSol.rho[0] = newSol.rho[1];
        newSol.rho[xN - 1] = bndCond.rhoRight;

        for (int i = 1; i < xN - 1; i++)
        {
            newSol.p[i] = pow((sol.rho[i + 1] + sol.rho[i - 1]) / 2, gamma);
        }
        newSol.p[0] = pow(sol.rho[0], gamma);
        newSol.p[xN - 1] = pow(sol.rho[xN - 1], gamma);

        sol = newSol;
    }
    return sol;
}

template<>
Solution Solver::Solve<LargeParticle>(const InitConditions& intiCond, const BndConditions& bndCond)
{
    Solution sol;
    sol.u = intiCond.u;
    sol.rho = intiCond.rho;
    sol.p = vector<double>(xN);
    for (int i = 0; i < xN; i++)
        sol.p[i] = pow(sol.rho[i], gamma);

    double s = tStep / xStep; 

    Smooth(sol, 0.1);

    for (int it = 0; it < tN; it++)
    {
        if (it % writingStep == 0)
        {
            cout << "t = " << it * tStep << "\n";
            double maxU = *max_element(sol.u.begin(), sol.u.end());
            cout << "CFL = " << tStep * maxU / xStep << "\n\n";

            WriteResults(to_string(it * tStep), sol);
        }

        Solution newSol;
        newSol.rho = vector<double>(xN);
        newSol.u = vector<double>(xN);
        newSol.p = vector<double>(xN);

        if (it % smoothingStep == 0)
            Smooth(sol, smoothingAlpha);

        // Eulerian step
        vector<double> uEul(xN);
        for (int i = 1; i < xN - 1; i++)
        {
            uEul[i] = sol.u[i] - s * 0.5 * (pow(sol.rho[i + 1], gamma) - pow(sol.rho[i - 1], gamma)) 
                / sol.rho[i];
        }
        uEul[0] = bndCond.uLeft;
        uEul[xN - 1] = bndCond.uRight;

        // Lagrangian step
        for (int i = 1; i < xN - 1; i++)
        {
            double rL = (i - 0.5) * xStep;
            double rR = (i + 0.5) * xStep;

            double deltaML = sol.rho[i - 1] * 0.5 * (uEul[i - 1] + uEul[i]) * tStep *
                4 * M_PI * pow(rL, 2);
            double deltaMR = sol.rho[i] * 0.5 * (uEul[i] + uEul[i + 1]) * tStep *
                4 * M_PI * pow(rR, 2);

            double dVol = 4.0 / 3 * M_PI * (pow(rR, 3) - pow(rL, 3));

            newSol.rho[i] = sol.rho[i] + (deltaML - deltaMR) / dVol;

            newSol.u[i] = (sol.rho[i] / newSol.rho[i]) * uEul[i] +
                (uEul[i - 1] * deltaML - uEul[i] * deltaMR) / (newSol.rho[i] * dVol);

            if (newSol.u[i] != newSol.u[i])
                throw runtime_error("nan appeared");
        }

        double deltaMR0 = sol.rho[0] * 0.5 * (uEul[0] + uEul[1]) * tStep *
            4 * M_PI * pow(0.5 * xStep, 2);
        double dV0 = 4.0 / 3 * M_PI * pow(0.5 * xStep, 3);
        newSol.rho[0] = sol.rho[0] - 2 * deltaMR0 / dV0;
        newSol.rho[xN - 1] = bndCond.rhoRight;

        newSol.u[0] = bndCond.uLeft;
        newSol.u[xN - 1] = bndCond.uRight;

        
        if (it % writingStep == 0)
        {
            ofstream out2;
            out2.open("delta_rho.txt");
            for (int i = 0; i < xN; i++)
                out2 << i * xStep << " " << newSol.rho[i] - sol.rho[i] << "\n";
            out2.close();
        }

        for (int i = 1; i < xN - 1; i++)
        {
            newSol.p[i] = pow((sol.rho[i + 1] + sol.rho[i - 1]) / 2, gamma);
        }
        newSol.p[0] = pow(sol.rho[0], gamma);
        newSol.p[xN - 1] = pow(sol.rho[xN - 1], gamma);

        sol = newSol;
    }
    return sol;
}