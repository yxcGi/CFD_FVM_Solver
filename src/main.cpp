#include <iostream>
#include "Laplacian.hpp"
#include "Vector.hpp"
#include "Geometry/Mesh.h"
#include "Tensor.hpp"
#include <Field.hpp>
#include "SparseMatrix.hpp"
#include <cmath>
#include "Solver.hpp"
#include "Div.hpp"
#include "Source.hpp"
#include "SIMPLE.hpp"
#include "SIMPLE1.hpp"





using Scalar = double;
int main() {

    // 三维方腔(非结构)
#if 1
     // 读取网格
    Mesh mesh("tempFile/OpenFOAM_tutorials/cavity3D_4/constant/polyMesh");

    // 定义速度场/压力场
    Field<Vector<Scalar>> U("U", &mesh);
    Field<Scalar> p("p", &mesh);

    // 初始化速度场/压力场
    U.setValue(Vector<Scalar>(0, 0, 0));
    p.setValue(0.0);

    // 设置边界条件
    U.setBoundaryCondition("movingWall",
                           1, 0, Vector<Scalar>(5, 0, 0));
    U.setBoundaryCondition("fixedWalls",
                           1, 0, Vector<Scalar>(0, 0, 0));

    p.setBoundaryCondition("movingWall", 0, 1, 0.0);
    p.setBoundaryCondition("fixedWalls", 0, 1, 0.0);


    // 设置密度/粘度
    FaceField<Scalar> rho("rho", &mesh);
    rho.setValue(1.0);

    FaceField<Scalar> nu("nu", &mesh);
    nu.setValue(0.01);

    // SIMPLE算法参数
    algorithm::simple::SIMPLE::Options options;
    options.maxOuterIterations = 5000;      // 外迭代次数
    options.alphaU = 0.7;                   // 松弛因子
    options.alphaP = 0.3;
    options.convergenceTolerance = 1e-8;
    options.nNonOrthogonalCorrectors = 2;
    options.divScheme = fvm::DivType::SUD;  // 目前只支持FUD
    options.useParallel = true;

    // 如果是压力出口，例如 patch 名为 "outlet"，则打开：
    // options.fixedPressurePatches = {"outlet"};

    algorithm::simple::SIMPLE solver(U, p, rho, nu, options);

    solver.solve();

    U.writeToFile("U_SIMPLE.dat");
    p.writeToFile("p_SIMPLE.dat");

    // std::cout << "Final continuity residual = " << residual.continuity << std::endl;

    return 0;
#endif


    // 台阶流
#if 0
    auto start = std::chrono::high_resolution_clock::now();
    // 读取网格
    Mesh mesh("tempFile/OpenFOAM_tutorials/pitzDailySteady/constant/polyMesh");

    // 定义速度场/压力场
    Field<Vector<Scalar>> U("U", &mesh);
    Field<Scalar> p("p", &mesh);

    // 初始化速度场/压力场
    U.setValue(Vector<Scalar>(0, 0, 0));
    p.setValue(0.0);

    // 设置边界条件
    U.setBoundaryCondition("inlet",
                           1, 0, Vector<Scalar>(10, 0, 0));
    U.setBoundaryCondition("outlet",
                           0, 1, Vector<Scalar>(0, 0, 0));
    U.setBoundaryCondition("upperWall",
                           1, 0, Vector<Scalar>(0, 0, 0));
    U.setBoundaryCondition("lowerWall",
                           1, 0, Vector<Scalar>(0, 0, 0));

    p.setBoundaryCondition("inlet",
                           0, 1, 0.0);
    p.setBoundaryCondition("outlet",
                           1, 0, 0.0);
    p.setBoundaryCondition("upperWall",
                           0, 1, 0.0);
    p.setBoundaryCondition("lowerWall",
                           0, 1, 0.0);

    // 设置密度/粘度
    FaceField<Scalar> rho("rho", &mesh);
    rho.setValue(1.0);

    FaceField<Scalar> nu("nu", &mesh);
    nu.setValue(0.01);

    // SIMPLE算法参数
    algorithm::simple::SIMPLE::Options options;
    options.maxOuterIterations = 500;      // 外迭代次数
    options.alphaU = 0.7;                   // 松弛因子
    options.alphaP = 0.3;
    options.convergenceTolerance = 1e-8;
    options.nNonOrthogonalCorrectors = 2;
    options.divScheme = fvm::DivType::SUD;  // 目前只支持FUD
    options.fixedPressurePatches = { "outlet" };
    options.useParallel = true;

    // 如果是压力出口，例如 patch 名为 "outlet"，则打开：
    // options.fixedPressurePatches = {"outlet"};

    algorithm::simple::SIMPLE solver(U, p, rho, nu, options);

    solver.solve();
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000 << "s" << std::endl;

    U.writeToFile("U_SIMPLE.dat");
    p.writeToFile("p_SIMPLE.dat");

    // std::cout << "Final continuity residual = " << residual.continuity << std::endl;

    return 0;
#endif


    // 二维方腔(非结构)
#if 0
    // 读取网格
    Mesh mesh("tempFile/OpenFOAM_tutorials/cavity2D_tri/constant/polyMesh");

    // 定义速度场/压力场
    Field<Vector<Scalar>> U("U", &mesh);
    Field<Scalar> p("p", &mesh);

    // 初始化速度场/压力场
    U.setValue(Vector<Scalar>(0, 0, 0));
    p.setValue(0.0);

    // 设置边界条件
    U.setBoundaryCondition("movingWall",
                           1, 0, Vector<Scalar>(5, 0, 0));
    U.setBoundaryCondition("fixedWalls",
                           1, 0, Vector<Scalar>(0, 0, 0));

    p.setBoundaryCondition("movingWall", 0, 1, 0.0);
    p.setBoundaryCondition("fixedWalls", 0, 1, 0.0);


    // 设置密度/粘度
    FaceField<Scalar> rho("rho", &mesh);
    rho.setValue(1.0);

    FaceField<Scalar> nu("nu", &mesh);
    nu.setValue(0.01);

    // SIMPLE算法参数
    algorithm::simple::SIMPLE::Options options;
    options.maxOuterIterations = 5000;      // 外迭代次数
    options.alphaU = 0.7;                   // 松弛因子
    options.alphaP = 0.3;
    options.convergenceTolerance = 1e-8;
    options.nNonOrthogonalCorrectors = 2;
    options.divScheme = fvm::DivType::SUD;  // 目前只支持FUD
    options.useParallel = true;

    // 如果是压力出口，例如 patch 名为 "outlet"，则打开：
    // options.fixedPressurePatches = {"outlet"};

    algorithm::simple::SIMPLE solver(U, p, rho, nu, options);

    solver.solve();

    U.writeToFile("U_SIMPLE.dat");
    p.writeToFile("p_SIMPLE.dat");

    // std::cout << "Final continuity residual = " << residual.continuity << std::endl;

    return 0;
#endif


    // 二维方腔(结构)
#if 0
    // 读取网格
    Mesh mesh("tempFile/OpenFOAM_tutorials/cavity/constant/polyMesh");

    // 定义速度场/压力场
    Field<Vector<Scalar>> U("U", &mesh);
    Field<Scalar> p("p", &mesh);

    // 初始化速度场/压力场
    U.setValue(Vector<Scalar>(0, 0, 0));
    p.setValue(0.0);

    // 设置边界条件
    U.setBoundaryCondition("movingWall",
                           1, 0, Vector<Scalar>(50, 0, 0));
    U.setBoundaryCondition("leftWalls",
                           1, 0, Vector<Scalar>(0, 0, 0));
    U.setBoundaryCondition("bottomWalls",
                           1, 0, Vector<Scalar>(0, 0, 0));
    U.setBoundaryCondition("rightWalls",
                           1, 0, Vector<Scalar>(0, 0, 0));

    p.setBoundaryCondition("movingWall",
                           0, 1, 0.0);
    p.setBoundaryCondition("leftWalls",
                           0, 1, 0.0);
    p.setBoundaryCondition("bottomWalls",
                           0, 1, 0.0);
    p.setBoundaryCondition("rightWalls",
                           0, 1, 0.0);

    // 设置密度/粘度
    FaceField<Scalar> rho("rho", &mesh);
    rho.setValue(1.0);

    FaceField<Scalar> nu("nu", &mesh);
    nu.setValue(0.01);

    // SIMPLE算法参数
    algorithm::simple::SIMPLE::Options options;
    options.maxOuterIterations = 5000;      // 外迭代次数
    options.alphaU = 0.7;                   // 松弛因子
    options.alphaP = 0.3;
    options.convergenceTolerance = 1e-8;
    options.nNonOrthogonalCorrectors = 2;
    options.divScheme = fvm::DivType::SUD;  // 目前只支持FUD
    options.useParallel = false;

    // 如果是压力出口，例如 patch 名为 "outlet"，则打开：
    // options.fixedPressurePatches = {"outlet"};

    algorithm::simple::SIMPLE solver(U, p, rho, nu, options);

    solver.solve();

    U.writeToFile("U_SIMPLE.dat");
    p.writeToFile("p_SIMPLE.dat");

    // std::cout << "Final continuity residual = " << residual.continuity << std::endl;

    return 0;
#endif




    // 源项测试
#if 0
    using Scalar = double;
    // 读取网格
    Mesh mesh("tempFile/OpenFOAM_tutorials/cavity/constant/polyMesh");

    // 创建标量场
    Field<Scalar> phi("T", &mesh);
    Field<Scalar> p("p", &mesh);
    phi.setValue(500);
    p.setValue(0);

    phi.setBoundaryCondition("movingWall", 0, 1, 0);
    phi.setBoundaryCondition("leftWalls", 1, 0, 100);
    phi.setBoundaryCondition("bottomWalls", 1, 0, 0);
    phi.setBoundaryCondition("rightWalls", 0, 1, 0);
    phi.cellToFace();       // 若是第一步，只是将边界面的场根据边界条件进行更新

    // 稀疏矩阵
    FaceField<Scalar> gamma("gamma", &mesh);
    gamma.setValue(20);
    SparseMatrix<Scalar> A_b(&mesh);

    fvm::Source(A_b, p.getCellField());
#endif

#if 0
    try {
        Mesh mesh("tempFile/OpenFOAM_tutorials/pitzDailySteady/constant/polyMesh");

        SparseMatrix<Scalar> A_b(&mesh);
        A_b.setValue(0, 0, 99);
        A_b.setValue(0, 1, 99);
        // A_b.setValue(0, 2, 99);
        // A_b.setValue(0, 3, 99);
        // A_b.setValue(0, 4, 99);

        // 打印第一行前5元素
        for (int i = 0; i < 5; i++) {
            std::cout << A_b.at(0, i) << " ";
        }
        std::cout << std::endl;

    }
    catch (std::exception& e) {
        std::cerr << "Exception: " << e.what() << std::endl;
    }
#endif


#if 0
    try {
        using Scalar = double;
        std::vector<std::vector<Scalar>> B{
            { 0, 0, 0, 0, 0 },
            { 0, 1, 0, 0, 0 },
            { 0, 0, 0, 0, 0 },
            { 2, 0, 0, 0, 0 },
            { 0, 0, 0, 4, 0 }
        };

        std::vector<std::vector<Scalar>> A{
            { 1, 1 },
            { 1, 2 }
        };
        std::vector<Scalar> b{ 3, 5 };


        SparseMatrix<Scalar> B_sparse(A);
        B_sparse.setB(b);
        B_sparse.printMatrix();
        B_sparse.compress();

        Solver<Scalar> solver(B_sparse, Solver<Scalar>::Method::Jacobi, 100);
        solver.init(b);
        solver.solve();

        // B_sparse.setValue(0, 0, 99);
        // B_sparse.setValue(0, 1, 99);
        // B_sparse.setValue(1, 0, 99);
        // B_sparse.setValue(1, 1, 99);

        B_sparse.printMatrix();


    }
    catch (std::exception& e) {
        std::cerr << "Exception: " << e.what() << std::endl;
    }
#endif


    // 对流项测试程序
#if 0
    try {
        using Scalar = double;
        // 读取网格
        Mesh mesh("tempFile/OpenFOAM_tutorials/cavity/constant/polyMesh");

        // 创建标量场
        Field<Scalar> phi("T", &mesh);

        phi.setValue(
            [](Scalar x, Scalar y, Scalar) {
                if (y > x) {
                    return 100;
                }
                else {
                    return 0;
                }
            }
        );

        phi.setBoundaryCondition("movingWall", 0, 1, 0);
        phi.setBoundaryCondition("leftWalls", 1, 0, 100);
        phi.setBoundaryCondition("bottomWalls", 1, 0, 0);
        phi.setBoundaryCondition("rightWalls", 0, 1, 0);
        phi.cellToFace();       // 若是第一步，只是将边界面的场根据边界条件进行更新

        // 稀疏矩阵
        FaceField<Scalar> rho("rho", &mesh);
        rho.setValue(1);


        FaceField<Vector<Scalar>> U("U", &mesh);
        U.setValue(
            [](Scalar x, Scalar, Scalar) {
                return Vector<Scalar>(10, 10, 0);
            }
        );


        SparseMatrix<Scalar> A_b(&mesh);

        // 创建稀疏矩阵
        // 对于非第一类边界条件需要循环迭代才可求解
        for (int i = 0; i < 1000; i++) {
            fvm::Div(A_b, rho, phi, U, fvm::DivType::FUD);
            // A_b.printMatrix();

            // getchar();

            Solver<Scalar> solver(A_b, Solver<Scalar>::Method::Jacobi, 100000);

            solver.init(phi.getCellField().getData());
            solver.setTolerance(1e-15);
            solver.relax(0.5);

            Scalar residual = solver.Error();
            std::cout << "residual: " << residual << " " << i << std::endl;
            if (residual < 1e-6) {
                break;
            }

            // solver.setParallel();

            // 计算求解时间
            // auto start = std::chrono::high_resolution_clock::now();
            solver.solve();
            // auto end = std::chrono::high_resolution_clock::now();
            // std::cout << "计算耗时：" << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;
            phi.cellToFace();
            A_b.clear();

        }

        phi.writeToFile("phi.dat");

        // getchar();
    }
    catch (const std::exception& e) {
        std::cerr << "Exception: " << e.what() << std::endl;
        return 1;
    }
#endif 


    // 矢量求解测试
#if 0
    try {
        using Scalar = double;
        // 读取网格
        Mesh mesh("tempFile/OpenFOAM_tutorials/cavity/constant/polyMesh");

        // 创建标量场
        Field<Vector<Scalar>> phi("U", &mesh);

        phi.setValue(Vector<Scalar>(0, 0, 0));

        // phi.setBoundaryCondition("movingWall", 1, 0, Vector<Scalar>(100, 0, 0));
        // phi.setBoundaryCondition("leftWalls", 1, 0, Vector<Scalar>(0, 0, 0));
        // phi.setBoundaryCondition("bottomWalls", 1, 0, Vector<Scalar>(0, 0, 0));
        // phi.setBoundaryCondition("rightWalls", 1, 0, Vector<Scalar>(0, 0, 0));
        phi.setBoundaryCondition("movingWall", 0, 1, Vector<Scalar>(0, 0, 0));
        phi.setBoundaryCondition("leftWalls", 1, 0, Vector<Scalar>(100, 0, 0));
        phi.setBoundaryCondition("bottomWalls", 1, 0, Vector<Scalar>(0, 0, 0));
        phi.setBoundaryCondition("rightWalls", 0, 1, Vector<Scalar>(0, 0, 0));
        phi.cellToFace();       // 若是第一步，只是将边界面的场根据边界条件进行更新

        // 稀疏矩阵
        SparseMatrix<Vector<Scalar>> A_b(&mesh);
        FaceField<Scalar> gamma("gamma", &mesh);
        gamma.setValue(20);

        // 创建稀疏矩阵
        // 对于非第一类边界条件需要循环迭代才可求解
        for (int i = 0; i < 1000; i++) {
            fvm::Laplacian(A_b, gamma, phi);

            Solver<Vector<Scalar>> solver(A_b, Solver<Vector<Scalar>>::Method::Jacobi, 100000);

            Scalar residual = solver.Error();
            std::cout << "residual: " << residual << " " << i << std::endl;
            if (residual < 1e-6) {
                break;
            }

            // solver.setParallel();

            // 计算求解时间
            // auto start = std::chrono::high_resolution_clock::now();
            solver.solve();
            // auto end = std::chrono::high_resolution_clock::now();
            // std::cout << "计算耗时：" << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;
            phi.cellToFace();
            A_b.clear();
        }
        phi.writeToFile("phi.dat");

        // getchar();
    }
    catch (const std::exception& e) {
        std::cerr << "Exception: " << e.what() << std::endl;
        return 1;
    }
#endif 


    // 方腔第二类边界条件测试
#if 0
    try {
        using Scalar = double;
        // 读取网格
        Mesh mesh("tempFile/OpenFOAM_tutorials/cavity/constant/polyMesh");

        // 创建标量场
        Field<Scalar> phi("T", &mesh);

        phi.setValue(500);

        phi.setBoundaryCondition("movingWall", 0, 1, 0);
        phi.setBoundaryCondition("leftWalls", 1, 0, 100);
        phi.setBoundaryCondition("bottomWalls", 1, 0, 0);
        phi.setBoundaryCondition("rightWalls", 0, 1, 0);
        phi.cellToFace();       // 若是第一步，只是将边界面的场根据边界条件进行更新

        // 稀疏矩阵
        FaceField<Scalar> gamma("gamma", &mesh);
        gamma.setValue(20);
        SparseMatrix<Scalar> A_b(&mesh);

        // 创建稀疏矩阵
        // 对于非第一类边界条件需要循环迭代才可求解
        for (int i = 0; i < 1; i++) {
            fvm::Laplacian(A_b, gamma, phi);

            Solver<Scalar> solver(A_b, Solver<Scalar>::Method::Jacobi, 100000);

            solver.init(phi.getCellField().getData());
            Scalar residual = solver.Error();
            std::cout << "residual: " << residual << " " << i << std::endl;
            if (residual < 1e-6) {
                break;
            }

            // solver.setParallel();

            // 计算求解时间
            // auto start = std::chrono::high_resolution_clock::now();
            solver.solve();
            // auto end = std::chrono::high_resolution_clock::now();
            // std::cout << "计算耗时：" << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;
            phi.cellToFace();
            A_b.clear();

        }

        phi.writeToFile("phi.dat");

        // getchar();
    }
    catch (const std::exception& e) {
        std::cerr << "Exception: " << e.what() << std::endl;
        return 1;
    }
#endif 


    // 方腔第一类边界条件测试
#if 0
    try {
        using Scalar = double;
        // 读取网格
        Mesh mesh("tempFile/OpenFOAM_tutorials/cavity/constant/polyMesh");

        // 创建标量场
        Field<Scalar> phi("T", &mesh);


        phi.setValue(500);

        phi.setBoundaryCondition("movingWall", 0, 1, 0);
        phi.setBoundaryCondition("leftWalls", 1, 0, 0);
        phi.setBoundaryCondition("bottomWalls", 1, 0, 0);
        phi.setBoundaryCondition("rightWalls", 1, 0, 0);
        phi.cellToFace();       // 若是第一步，只是将边界面的场根据边界条件进行更新

        // 稀疏矩阵
        FaceField<Scalar> gamma("gamma", &mesh);
        gamma.setValue(20);
        SparseMatrix<Scalar> A_b(&mesh);

        // 创建稀疏矩阵
        // 对于非第一类边界条件需要循环迭代才可求解
        for (int i = 0; i < 1000; i++) {
            fvm::Laplacian(A_b, gamma, phi);

            Solver<Scalar> solver(A_b, Solver<Scalar>::Method::Jacobi, 100000);

            solver.init(phi.getCellField_0().getData());
            Scalar residual = solver.Error();
            std::cout << "residual: " << residual << " " << i << std::endl;
            if (residual < 1e-6) {
                break;
            }

            // solver.setParallel();

            // 计算求解时间
            // auto start = std::chrono::high_resolution_clock::now();
            solver.solve();
            // auto end = std::chrono::high_resolution_clock::now();
            // std::cout << "计算耗时：" << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;
            phi.cellToFace();
            A_b.clear();

        }

        phi.writeToFile("phi.dat");

        // getchar();
    }
    catch (const std::exception& e) {
        std::cerr << "Exception: " << e.what() << std::endl;
        return 1;
    }
#endif 


    // 带源项的导热问题求解
#if 0
    try {
        using Scalar = double;
        // 读取网格
        Mesh mesh("tempFile/OpenFOAM_tutorials/cavity/constant/polyMesh");

        // 创建标量场
        Field<Scalar> phi("T", &mesh);
        CellField<Scalar> Q("Q", &mesh);        // 源项
        Q.setValue(
            [](Scalar x, Scalar y, Scalar z) {
                if ((x - 0.05) * (x - 0.05) + (y - 0.05) * (y - 0.05) < 0.0001) {
                    return 100000;
                }
                else {
                    return 0;
                }
            }
        );

        phi.setValue(500);

        phi.setBoundaryCondition("movingWall", 0, 1, 0);
        phi.setBoundaryCondition("leftWalls", 1, 0, 0);
        phi.setBoundaryCondition("bottomWalls", 1, 0, 0);
        phi.setBoundaryCondition("rightWalls", 1, 0, 0);
        phi.cellToFace();       // 若是第一步，只是将边界面的场根据边界条件进行更新

        // 稀疏矩阵
        FaceField<Scalar> gamma("gamma", &mesh);
        gamma.setValue(20);
        SparseMatrix<Scalar> A_b(&mesh);

        // 创建稀疏矩阵
        // 对于非第一类边界条件需要循环迭代才可求解
        for (int i = 0; i < 1000; i++) {
            fvm::Laplacian(A_b, gamma, phi);
            fvm::Source(A_b, Q);

            Solver<Scalar> solver(A_b, Solver<Scalar>::Method::Jacobi, 100000);

            solver.init(phi.getCellField_0().getData());
            Scalar residual = solver.Error();
            std::cout << "residual: " << residual << " " << i << std::endl;
            if (residual < 1e-6) {
                break;
            }

            // solver.setParallel();

            // 计算求解时间
            // auto start = std::chrono::high_resolution_clock::now();
            solver.solve();
            // auto end = std::chrono::high_resolution_clock::now();
            // std::cout << "计算耗时：" << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;
            phi.cellToFace();
            A_b.clear();

        }

        phi.writeToFile("phi.dat");

        // getchar();
    }
    catch (const std::exception& e) {
        std::cerr << "Exception: " << e.what() << std::endl;
        return 1;
    }
#endif 




    // cellField.cellToFace();
    // phi.setBoundaryCondition("const std::string &name", Scalar a, Scalar b, const double &c)


    // // 测试

    // using namespace std;

    // // 测试 Vector
    // Vector<double> v1(1.0, 2.0, 3.0);
    // Vector<double> v2(4.0, 5.0, 6.0);

    // cout << "v1 = " << v1 << endl;
    // cout << "v2 = " << v2 << endl;

    // cout << "v1 + v2 = " << v1 + v2 << endl;
    // cout << "v1 - v2 = " << v1 - v2 << endl;
    // cout << "v1 * 2 = " << v1 * 2.0 << endl;
    // cout << "v1 · v2 = " << (v1 & v2) << endl;
    // cout << "v1 × v2 = " << (v1 ^ v2) << endl;
    // cout << "|v1| = " << v1.magnitude() << endl;
    // cout << "unit(v1) = " << v1.unitVector() << endl;

    // // 测试 Tensor
    // Tensor<double> t1(1, 2, 3, 4, 5, 6, 7, 8, 9);
    // Tensor<double> t2(9, 8, 7, 6, 5, 4, 3, 2, 1);

    // cout << "t1 = " << t1 << endl;
    // cout << "t2 = " << t2 << endl;
    // cout << "t1 + t2 = " << t1 + t2 << endl;
    // cout << "t1 * v1 = " << t1 * v1 << endl;

    // // 测试并矢
    // Tensor<Scalar> dyad = v1 * v2;
    // cout << "v1 ⊗ v2 = " << dyad << endl;

    // return 0;
    // }
}