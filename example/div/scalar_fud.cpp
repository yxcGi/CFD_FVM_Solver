#include <exception>
#include <iostream>

#include "Div.hpp"
#include "Field.hpp"
#include "Geometry/Mesh.h"
#include "Solver.hpp"
#include "SparseMatrix.hpp"
#include "Vector.hpp"

using Scalar = double;

int main()
{
    try {
        Mesh mesh("tempFile/OpenFOAM_tutorials/cavity/constant/polyMesh");

        Field<Scalar> phi("T", &mesh);

        phi.setValue(
            [](Scalar x, Scalar y, Scalar) {
                return (y > x) ? 100 : 0;
            }
        );

        phi.setBoundaryCondition("movingWall", 0, 1, 0);
        phi.setBoundaryCondition("leftWalls", 1, 0, 100);
        phi.setBoundaryCondition("bottomWalls", 1, 0, 0);
        phi.setBoundaryCondition("rightWalls", 0, 1, 0);
        phi.cellToFace();

        FaceField<Scalar> rho("rho", &mesh);
        rho.setValue(1);

        FaceField<Vector<Scalar>> U("U", &mesh);
        U.setValue(
            [](Scalar, Scalar, Scalar) {
                return Vector<Scalar>(10, 10, 0);
            }
        );

        SparseMatrix<Scalar> A_b(&mesh);

        fvm::Div(A_b, rho, phi, U, fvm::DivType::FUD);

        Solver<Scalar> solver(A_b, Solver<Scalar>::Method::Jacobi, 100000);
        solver.init(phi.getCellField().getData());
        solver.setTolerance(1e-15);
        solver.relax(0.5);



        solver.solve();
        phi.cellToFace();

        phi.writeToFile("phi.dat");
    }
    catch (const std::exception& e) {
        std::cerr << "Exception: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
