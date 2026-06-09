#include <iostream>

#include "Field.hpp"
#include "Geometry/Mesh.h"
#include "SIMPLE.hpp"
#include "Vector.hpp"

using Scalar = double;

int main()
{
    // 三维方腔：非结构网格
    Mesh mesh("tempFile/OpenFOAM_tutorials/cavity3D_4/constant/polyMesh");

    Field<Vector<Scalar>> U("U", &mesh);
    Field<Scalar> p("p", &mesh);

    U.setValue(Vector<Scalar>(0, 0, 0));
    p.setValue(0.0);

    U.setBoundaryCondition("movingWall", 1, 0, Vector<Scalar>(5, 0, 0));
    U.setBoundaryCondition("fixedWalls", 1, 0, Vector<Scalar>(0, 0, 0));

    p.setBoundaryCondition("movingWall", 0, 1, 0.0);
    p.setBoundaryCondition("fixedWalls", 0, 1, 0.0);

    FaceField<Scalar> rho("rho", &mesh);
    rho.setValue(1.0);

    FaceField<Scalar> nu("nu", &mesh);
    nu.setValue(0.01);

    algorithm::simple::SIMPLE::Options options;
    options.maxOuterIterations = 50000;
    options.alphaU = 0.7;
    options.alphaP = 0.3;
    options.convergenceTolerance = 1e-8;
    options.nNonOrthogonalCorrectors = 2;
    options.divScheme = fvm::DivType::SUD;
    options.useParallel = true;

    algorithm::simple::SIMPLE solver(U, p, rho, nu, options);
    solver.solve();

    U.writeToFile("U_SIMPLE.dat");
    p.writeToFile("p_SIMPLE.dat");

    return 0;
}
