#include <exception>
#include <iostream>

#include "Field.hpp"
#include "Geometry/Mesh.h"
#include "Laplacian.hpp"
#include "Solver.hpp"
#include "SparseMatrix.hpp"
#include "Vector.hpp"

using Scalar = double;

int main()
{
    try {
        Mesh mesh("tempFile/OpenFOAM_tutorials/cavity/constant/polyMesh");

        Field<Vector<Scalar>> phi("U", &mesh);

        phi.setValue(Vector<Scalar>(0, 0, 0));

        phi.setBoundaryCondition("movingWall", 0, 1, Vector<Scalar>(0, 0, 0));
        phi.setBoundaryCondition("leftWalls", 1, 0, Vector<Scalar>(100, 0, 0));
        phi.setBoundaryCondition("bottomWalls", 1, 0, Vector<Scalar>(0, 0, 0));
        phi.setBoundaryCondition("rightWalls", 0, 1, Vector<Scalar>(0, 0, 0));
        phi.cellToFace();

        SparseMatrix<Vector<Scalar>> A_b(&mesh);

        FaceField<Scalar> gamma("gamma", &mesh);
        gamma.setValue(20);

        for (int i = 0; i < 1000; ++i) {
            fvm::Laplacian(A_b, gamma, phi);

            Solver<Vector<Scalar>> solver(
                A_b,
                Solver<Vector<Scalar>>::Method::Jacobi,
                100000
            );

            solver.init(phi.getCellField_0().getData());

            const Scalar residual = solver.Error();
            std::cout << "residual: " << residual << " " << i << std::endl;

            if (residual < 1e-6) {
                break;
            }

            solver.solve();
            phi.cellToFace();
            A_b.clear();
        }

        phi.writeToFile("phi.dat");
    }
    catch (const std::exception& e) {
        std::cerr << "Exception: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
