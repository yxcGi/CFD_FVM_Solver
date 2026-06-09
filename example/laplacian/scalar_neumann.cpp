#include <exception>
#include <iostream>

#include "Field.hpp"
#include "Geometry/Mesh.h"
#include "Laplacian.hpp"
#include "Solver.hpp"
#include "SparseMatrix.hpp"


int main() {
    using Scalar = double;
    Mesh mesh("tempFile/OpenFOAM_tutorials/cavity/constant/polyMesh");

    Field<Scalar> phi("T", &mesh);

    phi.setValue(500);

    phi.setBoundaryCondition("movingWall", 0, 1, 0);
    phi.setBoundaryCondition("leftWalls", 1, 0, 100);
    phi.setBoundaryCondition("bottomWalls", 1, 0, 0);
    phi.setBoundaryCondition("rightWalls", 0, 1, 0);
    phi.cellToFace();

    FaceField<Scalar> gamma("gamma", &mesh);
    gamma.setValue(20);

    SparseMatrix<Scalar> A_b(&mesh);

    fvm::Laplacian(A_b, gamma, phi);

    Solver<Scalar> solver(A_b, Solver<Scalar>::Method::Jacobi, 100000);
    solver.init(phi.getCellField().getData());


    solver.solve();
    phi.cellToFace();

    phi.writeToFile("phi.dat");

    return 0;
}
