#ifndef LAPLICIAN_H_
#define LAPLICIAN_H_

#include "SparseMatrix.hpp"
#include "Field.hpp"
#include <cassert>
#include "interpolation.hpp"

namespace fvm
{
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
            // assert(Ef_mag > 0);// 此处采用断言，Ef_mag一定大于0 测试用
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

        // 考虑边界面
        for (ULL faceId : boundaryFaceIndexes)
        {
            const Face& face = faces[faceId];   // 面
            Scalar faceArea = face.getArea();   // face面积
            const Vector<Scalar>& faceNormal = face.getNormal();    // 面法向量
            ULL owner = face.getOwnerIndex();
            const Vector<Scalar>& cb = face.getCenter() - cells[owner].getCenter();      // owner->neighbor向量
            const Vector<Scalar>& unit_cb = cb.unitVector();    // owner->neighbor单位向量

            Vector<Scalar> Sf = faceArea * faceNormal;  // face面向量

            // 计算系数
            Scalar gamma_f = gamma[faceId];
            Scalar d_cb = cb.magnitude();
            Scalar Ef_mag = faceArea / (faceNormal & unit_cb);
            // assert(Ef_mag > 0);// 此处采用断言，Ef_mag一定大于0 测试用
            Vector<Scalar> Ef = Ef_mag * cb.unitVector();
            Vector<Scalar> Tf = Sf - Ef;

            Scalar coefficient_E = gamma_f * Ef_mag / d_cb;
            Tp coefficient_T = gamma_f * faceGradientField[faceId] & Tf;

            // 添加矩阵元素
            matrix(owner, owner) += coefficient_E;
            matrix.addB(owner, coefficient_E * facefield[faceId] + coefficient_T);
        }
    }
}




#endif // LAPLICIAN_H_