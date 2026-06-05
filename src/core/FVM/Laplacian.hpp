#ifndef LAPLICIAN_H_
#define LAPLICIAN_H_

#include "SparseMatrix.hpp"
#include "Field.hpp"
#include <cassert>
#include "interpolation.hpp"

namespace fvm {
    /**
     * @brief 扩散项离散函数
     * @tparam Tp 场类型
     * @param matrix 待赋值矩阵
     * @param gamma 广义扩散系数
     * @param phi 待离散的场
     */
    template<typename Tp>
    void Laplacian(
        SparseMatrix<Tp>& matrix,   // 稀疏矩阵
        const FaceField<Scalar>& gamma,        // 扩散系数
        Field<Tp>& phi              // 场值，此处非常量引用，需要将地址传入matrix，后续会进行修改
    );


#pragma region 函数实现
    template<typename Tp>
    void Laplacian(
        SparseMatrix<Tp>& matrix,
        const FaceField<Scalar>& gamma,
        Field<Tp>& phi
    )
    {
        using GradType = decltype(Tp()* Vector<Scalar>());
        Mesh* mesh = matrix.getMesh();
        const CellField<GradType>& cellGradientField = phi.getCellGradientField();

        // 判断输入参数是否有效
        if (!matrix.isValid())
        {
            std::cerr << "Laplacian() Error: matrix is not valid." << std::endl;
            throw std::invalid_argument("matrix is not valid.");
        }
        if (!gamma.isValid())
        {
            std::cerr << "Laplacian() Error: gamma is not valid." << std::endl;
            throw std::invalid_argument("gamma is not valid.");
        }
        if (!phi.isValid())
        {
            std::cerr << "Laplacian() Error: phi is not valid." << std::endl;
            throw std::invalid_argument("phi is not valid.");
        }
        if (!cellGradientField.isValid())
        {
            std::cerr << "Laplacian() Error: gradientField is not valid." << std::endl;
            throw std::invalid_argument("gradientField is not valid.");
        }
        // 判断是否矩阵是否为空指针网格（稀疏矩阵只能由mesh构造）
        if (mesh == nullptr)
        {
            std::cerr << "Laplacian() Error: mesh is nullptr." << std::endl;
            throw std::invalid_argument("mesh is nullptr.");
        }
        // 判断输出参数是否为同一网格
        if (mesh != gamma.getMesh() ||
            mesh != phi.getMesh())
        {
            std::cerr << "Laplacian() Error: The mesh of the input fields parameters is not the same as the mesh of the matrix." << std::endl;
            throw std::invalid_argument("The mesh of the input fields parameters is not the same as the mesh of the matrix.");
        }

        // 给矩阵设置场，或者检查矩阵中的场是否与离散项一致
        if (matrix.fieldPtr_ == nullptr)
        {
            matrix.fieldPtr_ = &phi;
        }
        else    // 若已经被其他离散函数设置过
        {
            if (matrix.fieldPtr_ != &phi)
            {
                std::cerr << "Laplacian() Error: The field of the input matrix is not the same as the field of the input fields parameters." << std::endl;
                throw std::invalid_argument("The field of the input matrix is not the same as the field of the input fields parameters.");
            }
            // 否则矩阵中的场与输入场一致，什么都不发生
        }

        using ULL = unsigned long long;
        using LL = long long;
        const std::vector<Face>& faces = mesh->getFaces();
        const std::vector<Cell>& cells = mesh->getCells();
        const std::vector<ULL>& internalFaceIndexes = mesh->getInternalFaceIndexes();
        const std::vector<ULL>& boundaryFaceIndexes = mesh->getBoundaryFaceIndexes();
        const FaceField<Tp>& facefield = phi.getFaceField();
        util::Interpolation<GradType> interpolator; // 用于插值的函数对象

        // 获取单元梯度场
        // 计算内部面上的梯度场
        FaceField<GradType> faceGradientField("faceGradientField", mesh);
        faceGradientField.setValue(GradType{}); // 初始化使其有效
        for (ULL faceId : internalFaceIndexes)
        {
            const Face& face = faces[faceId];
            ULL owner = face.getOwnerIndex();
            LL neighbor = face.getNeighborIndex();
            const Point& faceCenter = face.getCenter();
            const Vector<Scalar>& faceNormal = face.getNormal();

            // 获取主单元邻单元的梯度值
            const GradType& ownerGrad = cellGradientField[owner];
            const GradType& neighborGrad = cellGradientField[neighbor];
            // 获取主单元与邻单元与面的距离
            Scalar ownerDistance = (faceCenter - cells[owner].getCenter()) & faceNormal;
            Scalar neighborDistance = (cells[neighbor].getCenter() - faceCenter) & faceNormal;
            Scalar alpha = ownerDistance / (ownerDistance + neighborDistance);

            // 插值
            faceGradientField[faceId] = interpolator(ownerGrad, neighborGrad, interpolation::Scheme::LINEAR, alpha);
        }
        // 边界面的梯度场直接用owner单元赋值
        for (ULL faceId : boundaryFaceIndexes)
        {
            faceGradientField[faceId] = cellGradientField[faces[faceId].getOwnerIndex()];
        }



        // 开始对面进行遍历，内部面与边界面分开，需考虑非正交修正，S_f = E_f + T_f
        // 先考虑内部面
        for (ULL faceId : internalFaceIndexes)
        {
            // 获取必要参数
            const Face& face = faces[faceId];   // 面
            Scalar faceArea = face.getArea();   // face面积
            const Vector<Scalar>& faceNormal = face.getNormal();    // 面法向量
            ULL owner = face.getOwnerIndex();
            LL neighbor = face.getNeighborIndex();
            const Vector<Scalar>& cb = cells[neighbor].getCenter() - cells[owner].getCenter();      // owner->neighbor向量
            const Vector<Scalar>& unit_cb = cb.unitVector();    // owner->neighbor单位向量

            Vector<Scalar> Sf = faceArea * faceNormal;  // face面向量


            // 计算系数
            Scalar gamma_f = gamma[faceId];
            Scalar d_cb = cb.magnitude();
            Scalar Ef_mag = faceArea / (faceNormal & unit_cb);
            Vector<Scalar> Ef = Ef_mag * cb.unitVector();
            Vector<Scalar> Tf = Sf - Ef;

            Scalar coefficient_E = gamma_f * Ef_mag / d_cb;
            Tp coefficient_T = gamma_f * faceGradientField[faceId] & Tf;

            // 添加矩阵元素
            matrix(owner, owner) += coefficient_E;
            matrix(owner, neighbor) -= coefficient_E;
            matrix(neighbor, neighbor) += coefficient_E;    // 这里是+号！！！！！
            matrix(neighbor, owner) -= coefficient_E;       // 这里是-号！！！！！


            // 考虑非正交修正，采用延迟修正
            matrix.addB(owner, coefficient_T);
            matrix.addB(neighbor, -coefficient_T);
        }

        // // 考虑边界面
        // for (ULL faceId : boundaryFaceIndexes)
        // {
        //     const Face& face = faces[faceId];   // 面
        //     Scalar faceArea = face.getArea();   // face面积
        //     const Vector<Scalar>& faceNormal = face.getNormal();    // 面法向量
        //     ULL owner = face.getOwnerIndex();
        //     const Vector<Scalar>& cb = face.getCenter() - cells[owner].getCenter();      // owner->neighbor向量
        //     const Vector<Scalar>& unit_cb = cb.unitVector();    // owner->neighbor单位向量

        //     Vector<Scalar> Sf = faceArea * faceNormal;  // face面向量

        //     // 计算系数
        //     Scalar gamma_f = gamma[faceId];
        //     Scalar d_cb = cb.magnitude();
        //     Scalar Ef_mag = faceArea / (faceNormal & unit_cb);
        //     Vector<Scalar> Ef = Ef_mag * cb.unitVector();
        //     Vector<Scalar> Tf = Sf - Ef;

        //     Scalar coefficient_E = gamma_f * Ef_mag / d_cb;
        //     Tp coefficient_T = gamma_f * faceGradientField[faceId] & Tf;

        //     // 添加矩阵元素
        //     matrix(owner, owner) += coefficient_E;
        //     matrix.addB(owner, coefficient_E * facefield[faceId] + coefficient_T);
        // }

        // 再遍历边界面，需考虑边界条件
        for (const auto& [name, bc] : phi.getBoundaryConditions())
        {
            if (bc.getType() ==
                BoundaryPatch::BoundaryType::EMPTY) // 不需要处理empty边界
            {
                continue;
            }

            ULL nFace = bc.getNFace();
            ULL startFace = bc.getStartFace();

            const ULL endFace = startFace + nFace;

            // for (ULL boundaryFaceIndex = startFace; boundaryFaceIndex < startFace + nFace; ++boundaryFaceIndex) {
            //     const Face& face = faces[boundaryFaceIndex]; Scalar faceArea = face.getArea();
            //     ULL ownerIndex = face.getOwnerIndex();
            //     const Cell& ownerCell = cells[ownerIndex];
            //     const Vector<Scalar>& faceCenter = face.getCenter();
            //     const Vector<Scalar>& ownerCenter = ownerCell.getCenter();
            //     const Vector<Scalar>& normal = face.getNormal();
            //     // 面法向量
            //     Vector<Scalar> V_CB = faceCenter - ownerCenter;
            //     // 计算中间量, normal = E + T 
            //     Vector<Scalar> E = V_CB / (V_CB & normal);
            //     Vector<Scalar> T = normal - E;
            //     Scalar E_magnitude = E.magnitude();
            //     Scalar distance_CB = ownerCenter.getDistance(faceCenter); Scalar a = bc.get_a();
            //     Scalar b = bc.get_b();
            //     const Tp& c = bc.get_c();
            //     // 计算c1, c2 
            //     Scalar c1 = (b * E_magnitude) / (a * distance_CB + b * E_magnitude);

            //     const GradType& cellGradient = cellGradientField[ownerIndex];
            //     Tp c2 = (c - (b * cellGradient & T)) * distance_CB / (a * distance_CB + b * E_magnitude);
            //     // 计算C1, C2 
            //     Scalar gamma_f = gamma[boundaryFaceIndex];
            //     Scalar EfMag = E_magnitude * faceArea;
            //     Scalar C1 = gamma_f * (1 - c1) * EfMag / distance_CB;
            //     Tp C2 = gamma_f * (c2 * EfMag / distance_CB + (cellGradient & (T * faceArea)));
            //     matrix(ownerIndex, ownerIndex) += C1;
            //     matrix.addB(ownerIndex, C2);
            // }

            for (ULL boundaryFaceIndex = startFace;
                 boundaryFaceIndex < endFace;
                 ++boundaryFaceIndex)
            {
                const Face& face = faces[boundaryFaceIndex];

                const ULL ownerIndex = face.getOwnerIndex();
                const Cell& ownerCell = cells[ownerIndex];

                const Vector<Scalar>& faceCenter = face.getCenter();
                const Vector<Scalar>& ownerCenter = ownerCell.getCenter();

                // 面面积 S_f
                const Scalar faceArea = face.getArea();

                // 单位面法向量 n_f
                // 注意：这里要求 normal 是单位法向量，而不是面积法向量。
                const Vector<Scalar>& normal = face.getNormal();

                // 从 owner cell center 指向 boundary face center 的向量 d = x_f - x_C
                const Vector<Scalar> d = faceCenter - ownerCenter;

                // d_n = d · n_f
                // 对正常边界面来说，d_n 应该大于 0。
                const Scalar dn = d & normal;


                // E = d / (d · n)
                // T = n - E
                // normal = E + T
                const Vector<Scalar> E = d / dn;
                const Vector<Scalar> T = normal - E;

                const Scalar a = bc.get_a();
                const Scalar b = bc.get_b();

                // 注意：这里要求 bc.get_c() 返回的是 const Tp&。
                // 如果 get_c() 返回临时对象，这里会产生悬空引用。
                const Tp& c = bc.get_c();

                const GradType& cellGradient = cellGradientField[ownerIndex];

                // 非正交修正项：grad(phi)_C · T
                //
                // 对 scalar transport：
                //     GradType = Vector<Scalar>
                //     Tp       = Scalar
                //
                // 对 vector transport：
                //     GradType = Tensor<Scalar>
                //     Tp       = Vector<Scalar>
                const Tp gradT = cellGradient & T;

                // Robin 边界条件分母：
                //
                //     a * phi_f + b * d(phi)/dn = c
                //
                // 离散化后自然出现 a * dn + b。
                const Scalar denom = a * dn + b;
                const Scalar gamma_f = gamma[boundaryFaceIndex];

                // C1 = gamma_f * a * S_f / (a * dn + b)
                const Scalar C1 = gamma_f * a * faceArea / denom;

                // C2 = gamma_f * S_f * (c + a * dn * gradT) / (a * dn + b)
                const Tp C2 = gamma_f * faceArea * (c + (a * dn) * gradT) / denom;

                matrix(ownerIndex, ownerIndex) += C1;
                matrix.addB(ownerIndex, C2);
            }
        }
    }
}




#endif // LAPLICIAN_H_