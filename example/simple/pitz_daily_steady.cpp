#include <chrono>
#include <iostream>

#include "Field.hpp"
#include "Geometry/Mesh.h"
#include "SIMPLE.hpp"
#include "Vector.hpp"

using Scalar = double;

int main()
{
    const auto start = std::chrono::high_resolution_clock::now();

    // 台阶流
    Mesh mesh("tempFile/OpenFOAM_tutorials/pitzDailySteady/constant/polyMesh");

    Field<Vector<Scalar>> U("U", &mesh);
    Field<Scalar> p("p", &mesh);

    U.setValue(Vector<Scalar>(0, 0, 0));
    p.setValue(0.0);

    U.setBoundaryCondition("inlet", 1, 0, Vector<Scalar>(10, 0, 0));
    U.setBoundaryCondition("outlet", 0, 1, Vector<Scalar>(0, 0, 0));
    U.setBoundaryCondition("upperWall", 1, 0, Vector<Scalar>(0, 0, 0));
    U.setBoundaryCondition("lowerWall", 1, 0, Vector<Scalar>(0, 0, 0));

    p.setBoundaryCondition("inlet", 0, 1, 0.0);
    p.setBoundaryCondition("outlet", 1, 0, 0.0);
    p.setBoundaryCondition("upperWall", 0, 1, 0.0);
    p.setBoundaryCondition("lowerWall", 0, 1, 0.0);

    FaceField<Scalar> rho("rho", &mesh);
    rho.setValue(1.0);

    FaceField<Scalar> nu("nu", &mesh);
    nu.setValue(0.01);

    algorithm::simple::SIMPLE::Options options;
    options.maxOuterIterations = 500;
    options.alphaU = 0.7;
    options.alphaP = 0.3;
    options.convergenceTolerance = 1e-8;
    options.nNonOrthogonalCorrectors = 2;
    options.divScheme = fvm::DivType::SUD;
    options.fixedPressurePatches = {"outlet"};
    options.useParallel = true;

    algorithm::simple::SIMPLE solver(U, p, rho, nu, options);
    solver.solve();

    const auto end = std::chrono::high_resolution_clock::now();
    std::cout
        << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0
        << "s"
        << std::endl;

    U.writeToFile("U_SIMPLE.dat");
    p.writeToFile("p_SIMPLE.dat");

    return 0;
}
