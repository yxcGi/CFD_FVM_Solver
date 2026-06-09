#include <iostream>

#include "Field.hpp"
#include "Geometry/Mesh.h"
#include "Source.hpp"
#include "SparseMatrix.hpp"


int main()
{
    using Scalar = double;
    // 源项矩阵装配测试
    Mesh mesh("tempFile/OpenFOAM_tutorials/cavity/constant/polyMesh");

    Field<Scalar> T("T", &mesh);
    Field<Scalar> p("p", &mesh);

    T.setValue(500);
    p.setValue(0);

    T.setBoundaryCondition("movingWall", 0, 1, 0);
    T.setBoundaryCondition("leftWalls", 1, 0, 100);
    T.setBoundaryCondition("bottomWalls", 1, 0, 0);
    T.setBoundaryCondition("rightWalls", 0, 1, 0);
    T.cellToFace();

    FaceField<Scalar> gamma("gamma", &mesh);
    gamma.setValue(20);
    SparseMatrix<Scalar> A_b(&mesh);
    fvm::Laplacian(A_b, gamma, T);

    fvm::Source(A_b, p.getCellField());

    std::cout << "Source assembly test finished." << std::endl;

    return 0;
}
