#include <exception>
#include <iostream>

#include "Geometry/Mesh.h"
#include "SparseMatrix.hpp"

using Scalar = double;

int main()
{
    try {
        Mesh mesh("tempFile/OpenFOAM_tutorials/pitzDailySteady/constant/polyMesh");

        SparseMatrix<Scalar> A_b(&mesh);
        A_b.setValue(0, 0, 99);
        A_b.setValue(0, 1, 99);

        for (int i = 0; i < 5; ++i) {
            std::cout << A_b.at(0, i) << " ";
        }
        std::cout << std::endl;
    }
    catch (const std::exception& e) {
        std::cerr << "Exception: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
