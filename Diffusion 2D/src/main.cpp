//-----------------------------------------------------------------------------------//
//      Soving the diffusion equation            u_t = d (u_xx + u_yy)               //
//      in the rectangle                         R = [0,L] x [0,L]                   //
//                                                                                   //
//      with Neumann boundary conditions         u_x(0,y,t) = l, u_x(L,y,t) = r      //
//                                               u_x(x,0,t) = b, u_x(x,L,t) = t      //
//                                               l,r,b,t = Const                     //
//                                                                                   //
//      and initial condition                    u(x,y,0) = u0(x,y)                  //                      
//-----------------------------------------------------------------------------------//

#include <iostream>
#include <fstream>
#include <vector>

#include "diff2d.h"

using namespace std;

void ExportSolution(UniformGrid grid, vector<double> u)
{
    ofstream outputGrid;
    outputGrid.open("grid_params.txt");
    outputGrid << grid.lx << " " << grid.ly << " " << grid.nx << " " << grid.ny;
    outputGrid.close();

    ofstream outputU;
    outputU.open("u.txt");
    for (int i = 0; i < grid.nx; i++)
    {
        for (int j = 0; j < grid.ny; j++)
        {
            outputU << u[i + j * grid.nx] << endl;
        }
    }
    outputU.close();
}

int main()
{
    double d = 0.25;

    double tEnd = 1.0;
    double tStep = 0.05;

    UniformGrid grid;
    grid.lx = 1.0;
    grid.ly = 1.0;
    grid.nx = 100;
    grid.ny = 100;

    NeumannBC bc;
    bc.left = -1.0;
    bc.bottom = 1.0;
    bc.right = 0.0;
    bc.top = 0.0;

    Diffusion solver(d, tStep, grid, bc);

    vector<double> u0(grid.nx * grid.ny, 0.0);
    //for (int i = 0; i < grid.nx; i++)
    //{
    //    for (int j = 0; j < grid.ny; j++)
    //    {
    //        if (i <= grid.nx / 2)
    //            u0[i + grid.nx * j] = 1;
    //    }
    //}

    vector<double> u = solver.Solve(u0, tEnd);

    ExportSolution(grid, u);
}