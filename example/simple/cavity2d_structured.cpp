#include <iostream>

#include "Field.hpp"
#include "Geometry/Mesh.h"
#include "SIMPLE.hpp"
#include "Vector.hpp"

using Scalar = double;

int main()
{
    // 二维方腔：结构网格
    Mesh mesh("tempFile/OpenFOAM_tutorials/cavity2D_4/constant/polyMesh");

    Field<Vector<Scalar>> U("U", &mesh);
    Field<Scalar> p("p", &mesh);

    U.setValue(Vector<Scalar>(0, 0, 0));
    p.setValue(0.0);

    // U.setBoundaryCondition("topWalls", 1, 0, Vector<Scalar>(4, 0, 0));  // Re = 400
    // U.setBoundaryCondition("topWalls", 1, 0, Vector<Scalar>(10, 0, 0));  // Re = 1000
    // U.setBoundaryCondition("topWalls", 1, 0, Vector<Scalar>(20, 0, 0));  // Re = 2000
    U.setBoundaryCondition("topWalls", 1, 0, Vector<Scalar>(50, 0, 0));  // Re = 5000

    U.setBoundaryCondition("leftWalls", 1, 0, Vector<Scalar>(0, 0, 0));
    U.setBoundaryCondition("bottomWalls", 1, 0, Vector<Scalar>(0, 0, 0));
    U.setBoundaryCondition("rightWalls", 1, 0, Vector<Scalar>(0, 0, 0));

    p.setBoundaryCondition("topWalls", 0, 1, 0.0);
    p.setBoundaryCondition("leftWalls", 0, 1, 0.0);
    p.setBoundaryCondition("bottomWalls", 0, 1, 0.0);
    p.setBoundaryCondition("rightWalls", 0, 1, 0.0);

    FaceField<Scalar> rho("rho", &mesh);
    rho.setValue(1.0);

    FaceField<Scalar> nu("nu", &mesh);
    nu.setValue(0.01);

    algorithm::simple::SIMPLE::Options options;
    options.maxOuterIterations = 50000;
    options.alphaU = 0.7;
    options.alphaP = 0.3;
    options.convergenceTolerance = 1e-8;
    options.nNonOrthogonalCorrectors = 3;
    options.divScheme = fvm::DivType::MUSCL;
    options.useParallel = true;

    algorithm::simple::SIMPLE solver(U, p, rho, nu, options);
    solver.solve();

    U.writeToFile("U_cavity2d_structured.dat");
    p.writeToFile("p_cavity2d_structured.dat");

    return 0;
}
