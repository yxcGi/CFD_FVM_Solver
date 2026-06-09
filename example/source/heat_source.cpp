#include <exception>
#include <iostream>

#include "Field.hpp"
#include "Geometry/Mesh.h"
#include "Laplacian.hpp"
#include "Solver.hpp"
#include "Source.hpp"
#include "SparseMatrix.hpp"

using Scalar = double;

int main()
{
    try {
        Mesh mesh("tempFile/OpenFOAM_tutorials/cavity/constant/polyMesh");

        Field<Scalar> phi("T", &mesh);


        phi.setValue(500);

        phi.setBoundaryCondition("movingWall", 0, 1, 0);
        phi.setBoundaryCondition("leftWalls", 1, 0, 0);
        phi.setBoundaryCondition("bottomWalls", 1, 0, 0);
        phi.setBoundaryCondition("rightWalls", 1, 0, 0);
        phi.cellToFace();

        FaceField<Scalar> gamma("gamma", &mesh);
        gamma.setValue(20);

        SparseMatrix<Scalar> A_b(&mesh);
        CellField<Scalar> Q("Q", &mesh);
        Q.setValue(10000);
        fvm::Laplacian(A_b, gamma, phi);
        fvm::Source(A_b, Q);

        Solver<Scalar> solver(A_b, Solver<Scalar>::Method::Jacobi, 100000);
        solver.init(phi.getCellField_0().getData());

        const Scalar residual = solver.Error();


        solver.solve();
        phi.cellToFace();
        A_b.clear();

        phi.writeToFile("phi.dat");
    }
    catch (const std::exception& e) {
        std::cerr << "Exception: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
