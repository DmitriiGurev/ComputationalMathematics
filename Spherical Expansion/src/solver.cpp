#include <solver.h>

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <array>
#include <cassert>

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

    out.open("u.txt");
    for (int i = 0; i < xN; i++)
        out << i * xStep << " " << sol.u[i] << "\n";
    out.close();
}

void Solver::WriteToVTK(string fileName, const Solution& sol)
{
    ofstream out;
    out.open(fileName + ".vtk");

    out << "# vtk DataFile Version 2.0\n";
    out << "Spherical Expansion\n";
    out << "ASCII\n";
    out << "DATASET UNSTRUCTURED_GRID\n";

    double L = 2;
    int n = 100;
    double h = 2 * L / (n - 1);

    int nTotal = n * n;

    out << "POINTS " << nTotal <<  " float\n";
    for (int i0 = 0; i0 < n; i0++)
    {
        for (int i1 = 0; i1 < n; i1++)
        {
                out << -L + h * i0 << " " <<
                       -L + h * i1 << " " << 
                       0.0 << "\n";
        }
    }
 
    out << "CELLS " << nTotal << " " << nTotal * 2 << "\n";
    for (int i0 = 0; i0 < n; i0++)
    {
        for (int i1 = 0; i1 < n; i1++)
        {
            out << 1 << " " << i1 + i0 * n << "\n";
        }
    } 

    out << "CELL_TYPES " << nTotal << "\n";
    for (int i = 0; i < nTotal; i++)
        out << 1 << "\n";
    
    out << "CELL_DATA " << nTotal << "\n";
    out << "SCALARS " << "density" << " double 1\n";
    out << "LOOKUP_TABLE default\n";
    for (int i0 = 0; i0 < n; i0++)
    {
        for (int i1 = 0; i1 < n; i1++)
        {
            double x0 = -L + h * i0;
            double x1 = -L + h * i1;

            double r = sqrt(x0 * x0 + x1 * x1);
            int ind = r / xStep;

            if (ind < sol.rho.size())
            {
                out << sol.rho[ind] << "\n";
            }
            else
            {
                out << 0.0 << "\n";
            }
        }
    }
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

struct Vec2
{
    array<double, 2> values;

    double operator[](int i) const
    {
        return values[i];
    }

    // Element-wise operations
    Vec2 operator+(const Vec2& v) const
    {
        return {values[0] + v.values[0], values[1] + v.values[1]};
    }

    Vec2 operator-(const Vec2& v) const
    {
        return {values[0] - v.values[0], values[1] - v.values[1]};
    }

    Vec2 operator*(const Vec2& v) const
    {
        return {values[0] * v.values[0], values[1] * v.values[1]};
    }

    Vec2 operator/(double d) const
    {
        return {values[0] / d, values[1] / d};
    }

    Vec2 operator*(double d) const
    {
        return {values[0] * d, values[1] * d};
    }

    Vec2 operator^(double d) const
    {
        return {pow(values[0], d), pow(values[1], d)};
    }

    friend ostream& operator<<(ostream& os, const Vec2& v);

    friend Vec2 operator/(double d, const Vec2& v);
    friend Vec2 operator/(const Vec2& n, const Vec2& d);

    friend Vec2 operator*(double d, const Vec2& v);
    friend Vec2 operator-(const Vec2& v);
};

ostream& operator<<(ostream& os, const Vec2& v)
{
    os << "Vector: {" << v.values[0] << ", " << v.values[1] << "}";
    return os;
}

Vec2 operator/(double d, const Vec2& v)
{
    return {d / v.values[0], d / v.values[1]};
}

Vec2 operator/(const Vec2& n, const Vec2& d)
{
    return {n.values[0] / d.values[0], n.values[1] / d.values[1]};
}

Vec2 operator*(double d, const Vec2& v)
{
    return v * d;
}

Vec2 operator-(const Vec2& v)
{
    return v * (-1);
}

template<>
Solution Solver::Solve<WenoLF>(const InitConditions& intiCond, const BndConditions& bndCond)
{
    Solution sol;
    sol.u = intiCond.u;
    sol.rho = intiCond.rho;
    sol.p = vector<double>(xN);
    for (int i = 0; i < xN; i++)
        sol.p[i] = pow(sol.rho[i], gamma);

    vector<Vec2> vC(xN);
    for (int i = 0; i < xN; i++) {
        double r = xStep * i;
        vC[i] = {pow(r, 2) * sol.rho[i], pow(r, 2) * sol.rho[i] * sol.u[i]};
    }

    auto f = [&](double r, const Vec2& v)
    {
        assert(v[0] > 0);
        return Vec2({
            v[1],
            pow(v[1], 2) / v[0] + pow(v[0], gamma) / pow(r, 2 * (gamma - 1))
            });
    };

    auto discretization = [&](const vector<Vec2>& vC)
    {
        // Flux splitting (Local LF)
        vector<Vec2> fPlus(xN);
        vector<Vec2> fMinus(xN);

        fPlus[0] = Vec2({0., 0.});
        fMinus[0] = Vec2({0., 0.});
        for (int i = 1; i < xN; i++) 
        {
            double r = xStep * i;

            // Maximum eigenvalue
            double s = max({
                abs(vC[i][1] / vC[i][0] +
                sqrt(gamma * pow(vC[i][0], gamma - 1) / pow(r, 2 * (gamma - 1)))),
                abs(vC[i][1] / vC[i][0] -
                sqrt(gamma * pow(vC[i][0], gamma - 1) / pow(r, 2 * (gamma - 1))))
                });

            fPlus[i] = 0.5 * (f(r, vC[i]) + s * vC[i]);
            fMinus[i] = 0.5 * (f(r, vC[i]) - s * vC[i]);
        }

        // Flux reconstruction      
        auto fP = [&](int ind)
        {
            if (ind < 0)
            {
                return fPlus[-ind];
            }
            else if (ind > xN - 1)
            {
                return fPlus[xN - 1];
            }
            else
            {
                return fPlus[ind];
            }
        };

        auto fM = [&](int ind)
        {
            if (ind < 0)
            {
                return fMinus[-ind];
            }
            else if (ind > xN - 1)
            {
                return fMinus[xN - 1];
            }
            else
            {
                return fMinus[ind];
            }
        };

        vector<Vec2> fBndR(xN - 1);
        vector<Vec2> fBndL(xN - 1);
        for (int i = 0; i < xN - 1; i++)
        {
            array<Vec2, 3> beta;
            array<Vec2, 3> alpha;
            array<Vec2, 3> weight;
            Vec2 eps = {1.0e-10, 1.0e-10};

            // R
            beta[0] = 13. / 12 * ((fP(i - 2) - 2 * fP(i - 1) + fP(i)) ^ 2.) +
                1. / 4 * ((fP(i - 2) - 4 * fP(i - 1) + 3 * fP(i)) ^ 2.);
            beta[1] = 13. / 12 * ((fP(i - 1) - 2 * fP(i) + fP(i + 1)) ^ 2.) +
                1. / 4 * ((fP(i - 1) - fP(i + 1)) ^ 2.);
            beta[2] = 13. / 12 * ((fP(i) - 2 * fP(i + 1) + fP(i + 2)) ^ 2.) +
                1. / 4 * ((3 * fP(i) - 4 * fP(i + 1) + fP(i + 2)) ^ 2.);

            alpha[0] = (1. / 10) / ((eps + beta[0]) ^ 2.);
            alpha[1] = (6. / 10) / ((eps + beta[1]) ^ 2.);
            alpha[2] = (3. / 10) / ((eps + beta[2]) ^ 2.);

            weight[0] = alpha[0] / (alpha[0] + alpha[1] + alpha[2]);
            weight[1] = alpha[1] / (alpha[0] + alpha[1] + alpha[2]);
            weight[2] = alpha[2] / (alpha[0] + alpha[1] + alpha[2]);

            fBndR[i] = 1. / 6 * weight[0] * (2 * fP(i - 2) - 7 * fP(i - 1) + 11 * fP(i)) +
                1. / 6 * weight[1] * (-fP(i - 1) + 5 * fP(i) + 2 * fP(i + 1)) +
                1. / 6 * weight[2] * (2 * fP(i) + 5 * fP(i + 1) - fP(i + 2));

            // L
            beta[0] = 13. / 12 * ((fM(i + 1) - 2 * fM(i + 2) + fM(i + 3)) ^ 2.) +
                1. / 4 * ((3 * fM(i + 1) - 4 * fM(i + 2) + fM(i + 3)) ^ 2.);
            beta[1] = 13. / 12 * ((fM(i) - 2 * fM(i + 1) + fM(i + 2)) ^ 2.) +
                1. / 4 * ((fM(i) - fM(i + 2)) ^ 2.);
            beta[2] = 13. / 12 * ((fM(i - 1) - 2 * fM(i) + fM(i + 1)) ^ 2.) +
                1. / 4 * ((fM(i - 1) - 4 * fM(i) + 3 * fM(i + 1)) ^ 2.);

            alpha[0] = (1. / 10) / ((eps + beta[0]) ^ 2.);
            alpha[1] = (6. / 10) / ((eps + beta[1]) ^ 2.);
            alpha[2] = (3. / 10) / ((eps + beta[2]) ^ 2.);

            weight[0] = alpha[0] / (alpha[0] + alpha[1] + alpha[2]);
            weight[1] = alpha[1] / (alpha[0] + alpha[1] + alpha[2]);
            weight[2] = alpha[2] / (alpha[0] + alpha[1] + alpha[2]);

            fBndL[i] = 1. / 6 * weight[2] * (-fM(i - 1) + 5 * fM(i) + 2 * fM(i + 1)) +
                1. / 6 * weight[1] * (2 * fM(i) + 5 * fM(i + 1) - fM(i + 2)) +
                1. / 6 * weight[0] * (11 * fM(i + 1) - 7 * fM(i + 2) + 2 * fM(i + 3));
        }

        vector<Vec2> rhs(xN);
        rhs[0] = {0., 0.};
        for (int i = 1; i < xN; i++)
            rhs[i] = {0., 2 * (i * xStep) * pow(vC[i][0] / pow(i * xStep, 2), gamma)};

        vector<Vec2> l(xN, Vec2({0., 0.}));
        for (int i = 1; i < xN - 1; i++)
            l[i] = rhs[i] - 1 / xStep * ((fBndR[i] - fBndR[i - 1]) + (fBndL[i] - fBndL[i - 1]));

        return l;
    };

    for (int it = 0; it < tN; it++)
    {
        if (it % 100 == 0)
            cout << it << "\n";

        if (it % writingStep == 0)
        {
            cout << "t = " << it * tStep << "\n";
            double maxU = *max_element(sol.u.begin(), sol.u.end());
            cout << "CFL = " << tStep * maxU / xStep << "\n\n";
            // WriteResults(to_string(it * tStep), sol);
            // WriteToVTK("vtk/sol_" + to_string(it * tStep), sol);
        }

        // Integration over time (RG-4)
        // (1)
        vector<Vec2> vC1 = vC;
        vector<Vec2> discr1 = discretization(vC);
        for (int i = 1; i < xN - 1; i++)
            vC1[i] = vC1[i] + 0.5 * tStep * discr1[i];
        
        // (2)
        vector<Vec2> vC2 = vC;
        vector<Vec2> discr2 = discretization(vC1);
        for (int i = 1; i < xN - 1; i++)
            vC2[i] = vC2[i] + 0.5 * tStep * discr2[i];

        // (3)
        vector<Vec2> vC3 = vC;
        vector<Vec2> discr3 = discretization(vC2);
        for (int i = 1; i < xN - 1; i++)
            vC3[i] = vC3[i] + tStep * discr3[i];

        // (4)
        vector<Vec2> newVC = vC;
        newVC[0] = {0., 0.};

        vector<Vec2> discr4 = discretization(vC3);
        for (int i = 1; i < xN - 1; i++)
            newVC[i] = newVC[i] +
                1. / 6 * tStep * (discr1[i] + 2 * discr2[i] + 2 * discr3[i] + 2 * discr4[i]);

        // TODO: Update BC
        newVC[xN - 1] = {
            pow(xStep * (xN - 1), 2) * bndCond.rhoRight,
            pow(xStep * (xN - 1), 2) * bndCond.rhoRight * bndCond.uRight
            };

        // Writing the solution
        if (it % writingStep == 0)
        {
            ofstream out;
            out.open("v.txt");
            if (out.is_open()) {
                for (int i = 0; i < xN; i++)
                    out << i * xStep << " " << vC[i][0] << " " << vC[i][1] << "\n";
            }
            out.close();
        }

        vC = newVC;

        Solution newSol;
        newSol.rho = vector<double>(xN);
        newSol.u = vector<double>(xN);
        newSol.p = vector<double>(xN);

        for (int i = 1; i < xN - 1; i++)
            newSol.rho[i] = newVC[i][0] / pow(i * xStep, 2);
        newSol.rho[0] = newSol.rho[1];
        newSol.rho[xN - 1] = bndCond.rhoRight;

        for (int i = 1; i < xN - 1; i++)
            newSol.u[i] = newVC[i][1] / newVC[i][0];
        newSol.u[0] = bndCond.uLeft;
        newSol.u[xN - 1] = bndCond.uRight;
        
        for (int i = 0; i < xN; i++)
            newSol.p[i] = pow(newSol.rho[i], gamma);

        sol = newSol;
    }

    return sol;
}