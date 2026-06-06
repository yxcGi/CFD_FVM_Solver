#include <iostream>

#include "Field.hpp"
#include "Geometry/Mesh.h"
#include "Source.hpp"
#include "SparseMatrix.hpp"

using Scalar = double;

int main()
{
    // 源项矩阵装配测试
    Mesh mesh("tempFile/OpenFOAM_tutorials/cavity/constant/polyMesh");

    Field<Scalar> phi("T", &mesh);
    Field<Scalar> p("p", &mesh);

    phi.setValue(500);
    p.setValue(0);

    phi.setBoundaryCondition("movingWall", 0, 1, 0);
    phi.setBoundaryCondition("leftWalls", 1, 0, 100);
    phi.setBoundaryCondition("bottomWalls", 1, 0, 0);
    phi.setBoundaryCondition("rightWalls", 0, 1, 0);
    phi.cellToFace();

    FaceField<Scalar> gamma("gamma", &mesh);
    gamma.setValue(20);

    SparseMatrix<Scalar> A_b(&mesh);

    fvm::Source(A_b, p.getCellField());

    std::cout << "Source assembly test finished." << std::endl;

    return 0;
}
