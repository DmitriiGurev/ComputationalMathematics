#include "tucker.h"
#include <iostream>

#include <eigen-3.4.0/unsupported/Eigen/CXX11/Tensor>
#include <eigen-3.4.0/Eigen/Dense>

using namespace Eigen;
using std::cout;

#include <iostream>

Tensor<double, 3> Folding3D(MatrixXd unfolding, int index, int I0, int I1, int I2)
{
    Tensor<double, 3> folding(I0, I1, I2);

    switch (index)
    {
    case 0:
        for (int i0 = 0; i0 < I0; i0++)
        {
            for (int i1 = 0; i1 < I1; i1++)
            {
                for (int i2 = 0; i2 < I2; i2++)
                {
                    folding(i0, i1, i2) = unfolding(i0, i2 + i1 * I2);
                }
            }
        }
        break;
    case 1:
        for (int i1 = 0; i1 < I1; i1++)
        {
            for (int i2 = 0; i2 < I2; i2++)
            {
                for (int i0 = 0; i0 < I0; i0++)
                {
                    folding(i0, i1, i2) = unfolding(i1, i0 + i2 * I0);
                }
            }
        }
        break;
    case 2:
        for (int i2 = 0; i2 < I2; i2++)
        {
            for (int i0 = 0; i0 < I0; i0++)
            {
                for (int i1 = 0; i1 < I1; i1++)
                {
                    folding(i0, i1, i2) = unfolding(i2, i1 + i0 * I1);
                }
            }
        }
        break;
    }
    return folding;
}

int main()
{
    Matrix<double, 3, 9, RowMajor> A0;
    A0 << 0.9073, 0.7158, -0.3698, 1.7842, 1.6970, 0.0151, 2.1236, -0.0740, 1.4429,
      0.8924, -0.4898, 2.4288, 1.7753, -1.5077, 4.0337, -0.6631, 1.9103, -1.7495,
      2.1488, 0.3054, 2.3753, 4.2495, 0.3207, 4.7146, 1.8260, 2.1335, -0.2716;

    Tensor<double, 3> tensor = Folding3D(A0, 0, 3, 3, 3);

    Tucker tucker(tensor);

    cout << "Accurate representation of the original tensor (eps = 0):\n" << tucker << "\n\n";

    cout << "Dimensions: " 
        << "{" << tucker.Dimensions()[0] << ", " << tucker.Dimensions()[1] << ", " << tucker.Dimensions()[2] << "}\n\n";

    cout << "Ranks: "
        << "{" << tucker.Ranks()[0] << ", " << tucker.Ranks()[1] << ", " << tucker.Ranks()[2] << "}\n\n";

    cout << "Factors:\n" 
        << tucker.U()[0] << "\n\n" << tucker.U()[1] << "\n\n" << tucker.U()[2] << "\n\n";

    cout << "Core:\n"
        << tucker.Core() << "\n\n";

    cout << "Frobenius norm:\n"
        << tucker.Norm() << "\n\n";

    cout << "--------------------------------------------------------------------------------------------------------------------\n\n";

    Tucker compressed(tensor, 0.1);
    cout << "Compressed version of the original tensor (eps = 0.1):\n" << compressed << "\n\n";

    cout << "Dimensions: "
        << "{" << compressed.Dimensions()[0] << ", " << compressed.Dimensions()[1] << ", " << compressed.Dimensions()[2] << "}\n\n";

    cout << "Ranks: "
        << "{" << compressed.Ranks()[0] << ", " << compressed.Ranks()[1] << ", " << compressed.Ranks()[2] << "}\n\n";

    cout << "Factors:\n"
        << compressed.U()[0] << "\n\n" << compressed.U()[1] << "\n\n" << compressed.U()[2] << "\n\n";

    cout << "Core:\n"
        << compressed.Core() << "\n\n";

    cout << "Frobenius norm:\n"
        << compressed.Norm() << "\n\n";
}
