#ifndef DIV_H_
#define DIV_H_

#include "Field.hpp"
#include "SparseMatrix.hpp"
#include "DivType.h"
#include "interpolation.hpp"
#include <cassert>


namespace fvm
{

    /**
     * @brief 对流项离散函数
     * @tparam Tp 场类型
     * @param matrix 待赋值矩阵
     * @param rho 密度
     * @param phi 待离散场
     * @param U 面速度
     * @param type 对流离散格式
     */
    template<typename Tp>
    void Div(
        SparseMatrix<Tp>& matrix,
        const FaceField<Scalar>& rho,   // 密度场
        Field<Tp>& phi,                 // 待离散场
        FaceField<Vector<Scalar>>& U,   // 速度场
        DivType type  // 离散格式
    );


    template<typename Tp>
    void Div(
        SparseMatrix<Tp>& matrix,
        const FaceField<Scalar>& rho,
        Field<Tp>& phi,
        FaceField<Vector<Scalar>>& U,
        DivType type)
    {
        Mesh* mesh = matrix.getMesh();

        // 检查输入参数是否有效
        if (!matrix.isValid())
        {
            std::cerr << "Div() Error: matrix is not valid." << std::endl;
            throw std::invalid_argument("matrix is not valid.");
        }
        if (!rho.isValid())
        {
            std::cerr << "Div() Error: rho is not valid." << std::endl;
            throw std::invalid_argument("rho is not valid.");
        }
        if (!phi.isValid())
        {
            std::cerr << "Div() Error: phi is not valid." << std::endl;
            throw std::invalid_argument("phi is not valid.");
        }

        // 判断矩阵是否为空指针网格（稀疏矩阵只能由mesh构造）
        if (mesh == nullptr)
        {
            std::cerr << "Div() Error: mesh is nullptr." << std::endl;
            throw std::invalid_argument("mesh is nullptr.");
        }
        // 参数是否为同一网格
        if (mesh != rho.getMesh() ||
            mesh != phi.getMesh() ||
            mesh != U.getMesh())
        {
            std::cerr << "Div() Error: The mesh of the input fields parameters is not the same as the mesh of the matrix." << std::endl;
            throw std::invalid_argument("The mesh of the input fields parameters is not the same as the mesh of the matrix.");
        }

        // 给矩阵设置场，或者检查矩阵中的场是否与离散项一致
        if (matrix.fieldPtr_ == nullptr)
        {
            matrix.fieldPtr_ = &phi;
        }
        else    // 若已被其他离散函数设置过
        {
            if (matrix.fieldPtr_ != &phi)
            {
                std::cerr << "Div() Error: The field of the matrix is not the same as the field of the input parameters." << std::endl;
                throw std::invalid_argument("The field of the matrix is not the same as the field of the input parameters.");
            }
            // 否则矩阵中的场与输入场一致，什么都不发生
        }

        // 获取必要参数
        using ULL = unsigned long long;
        using LL = long long;
        using GradType = decltype(Tp()* Vector<Scalar>());
        const std::vector<Face>& faces = mesh->getFaces();
        const std::vector<Cell>& cells = mesh->getCells();
        const std::vector<ULL>& internalFaceIndexes = mesh->getInternalFaceIndexes();
        const std::vector<ULL>& boundaryFaceIndexes = mesh->getBoundaryFaceIndexes();
        const FaceField<Tp>& phiFaceField = phi.getFaceField(); // 待离散场的面场
        const CellField<Tp>& phiCellField_0 = phi.getCellField();   // 获取本步场值
        const CellField<GradType>& cellGradient = phi.getCellGradientField();


        // 先计算每个面的质量流量mf，之后phi要左乘mf
        std::vector<Scalar> mf(mesh->getFaceNumber());

        for (ULL faceId = 0;
            faceId < mesh->getFaceNumber();
            ++faceId)
        {
            const Face& face = faces[faceId];
            const Vector<Scalar>& faceNormal = face.getNormal();

            // mf = rho * U · Sf
            mf[faceId] = rho[faceId] * (U[faceId] & face.getNormal() * face.getArea());
        }


        // 根据离散格式计算矩阵
        if (type == DivType::FUD)   // 一阶迎风
        {
            // 先遍历内部面
            for (ULL faceId : internalFaceIndexes)
            {
                const Face& face = faces[faceId];
                ULL owner = face.getOwnerIndex();
                LL neighbor = face.getNeighborIndex();
                Scalar m_f = mf[faceId];

                if (m_f >= 0)
                {
                    matrix(owner, owner) += m_f;
                    matrix(neighbor, owner) -= m_f;
                }
                else
                {
                    matrix(owner, neighbor) += m_f;
                    matrix(neighbor, neighbor) -= m_f;
                }
            }
            // 再遍历边界面
            for (ULL faceId : boundaryFaceIndexes)
            {
                const Face& face = faces[faceId];
                ULL owner = face.getOwnerIndex();
                Scalar m_f = mf[faceId];

                if (m_f >= 0)   // 用边界单元的
                {
                    matrix(owner, owner) += m_f;
                }
                else    // 如果流量从边界流入，直接用边界面的
                {
                    matrix.addB(owner, -phiFaceField[faceId] * m_f);
                }
            }
        }
        else if (type == DivType::SUD)  // 二阶迎风：一阶迎风+二阶延迟修正
        {
            // 先遍历内部面
            for (ULL faceId : internalFaceIndexes)
            {
                const Face& face = faces[faceId];
                ULL owner = face.getOwnerIndex();
                LL neighbor = face.getNeighborIndex();
                Scalar m_f = mf[faceId];
                Vector<Scalar> CD = cells[neighbor].getCenter() - cells[owner].getCenter();
                GradType ownerGrad = cellGradient[owner];
                GradType neighborGrad = cellGradient[neighbor];

                if (m_f >= 0)
                {
                    // matrix(owner, owner) += m_f * 1.5;
                    // matrix(owner, neighbor) -= m_f * 0.5;
                    // matrix.addB(owner, -m_f * (ownerGrad & CD));

                    // matrix(neighbor, neighbor) += m_f * 0.5;
                    // matrix(neighbor, owner) -= m_f * 1.5;
                    // matrix.addB(neighbor, m_f * (ownerGrad & CD));

                    // 一阶部分
                    matrix(owner, owner) += m_f;
                    matrix(neighbor, owner) -= m_f;

                    // 二阶修正部分 (二阶 - 一阶)
                    matrix.addB(owner, -m_f * 1.5 * phiCellField_0[owner] + m_f * 0.5 * phiCellField_0[neighbor] - m_f * (ownerGrad & CD) + m_f * phiCellField_0[owner]);

                    matrix.addB(neighbor, -m_f * 0.5 * phiCellField_0[neighbor] + m_f * 1.5 * phiCellField_0[owner] + m_f * (ownerGrad & CD) - m_f * phiCellField_0[owner]);
                }
                else
                {
                    // matrix(owner, owner) -= m_f * 0.5;
                    // matrix(owner, neighbor) += m_f * 1.5;
                    // matrix.addB(owner, m_f * (neighborGrad & CD));  // 这里DC = -CD，省略了

                    // matrix(neighbor, neighbor) += m_f * 1.5;
                    // matrix(neighbor, owner) -= m_f * 0.5;
                    // matrix.addB(neighbor, -m_f * (neighborGrad & CD));  // 这里DC = -CD，省略了

                    // 一阶部分
                    matrix(owner, neighbor) += m_f;
                    matrix(neighbor, neighbor) -= m_f;

                    // 二阶修正部分
                    matrix.addB(owner, +m_f * 0.5 * phiCellField_0[owner] - m_f * 1.5 * phiCellField_0[neighbor] + m_f * (neighborGrad & CD) + m_f * phiCellField_0[neighbor]);
                    matrix.addB(neighbor, +m_f * 1.5 * phiCellField_0[neighbor] - m_f * 0.5 * phiCellField_0[owner] - m_f * (neighborGrad & CD) - m_f * phiCellField_0[neighbor]);

                }
            }

            // 再遍历边界面
            for (ULL faceId : boundaryFaceIndexes)
            {
                const Face& face = faces[faceId];
                ULL owner = face.getOwnerIndex();
                Scalar m_f = mf[faceId];
                Vector<Scalar> Cf = face.getCenter() - cells[owner].getCenter();
                GradType ownerGrad = cellGradient[owner];

                // assert(Cf & face.getNormal() > 0);

                if (m_f >= 0)
                {
                    // matrix(owner, owner) += m_f * 1.5;
                    // matrix.addB(owner, m_f * (0.5 * phiFaceField[faceId] - 1.5 * (ownerGrad & Cf)));

                    // 一阶部分
                    matrix(owner, owner) += m_f;

                    // 二阶部分
                    // matrix.addB(owner, -m_f * 1.5 * phiCellField_0[owner] + m_f * (0.5 * phiFaceField[faceId] - 1.5 * (ownerGrad & Cf)) + m_f * phiCellField_0[owner]);
                }
                else
                {
                    matrix.addB(owner, -m_f * phiFaceField[faceId]);
                }
            }
        }
        else if (type == DivType::CD)   // 中心差分,等距结构化网格不可使用
        {
            const std::vector<Cell>& cells = mesh->getCells();
            // 内部面
            for (ULL faceId : internalFaceIndexes)  // 内部面
            {
                const Face& face = faces[faceId];
                ULL owner = face.getOwnerIndex();
                LL neighbor = face.getNeighborIndex();
                Scalar m_f = mf[faceId];
                Point faceCenter = face.getCenter();

                std::cout << "m_f: " << m_f << "\n";

                // assert(std::abs(m_f) > 1e-15); // 中心差分时，mf不应为0 测试用

                Scalar ownerDistance = (faceCenter - cells[owner].getCenter()).magnitude();
                Scalar neighborDistance = (cells[neighbor].getCenter() - faceCenter).magnitude();


                // 计算系数
                Scalar ownerCoefficient = m_f * neighborDistance / (ownerDistance + neighborDistance);
                Scalar neighborCoefficient = m_f * ownerDistance / (ownerDistance + neighborDistance);
                matrix(owner, owner) += ownerCoefficient;
                matrix(owner, neighbor) += neighborCoefficient;
                matrix(neighbor, owner) -= ownerCoefficient;
                matrix(neighbor, neighbor) -= neighborCoefficient;
            }

            // // 测试用，检查对角线元素是否为0
            // for (int i = 0; i < mesh->getCellNumber(); ++i)
            // {
            //     if (std::abs(matrix(i, i)) < 1e-15)
            //     {
            //         matrix(i, i) = 1e-15; // 防止对角线为0
            //         std::cout << "Div() Warning: The diagonal element of matrix is too small." << std::endl;
            //     }
            // }

            // 边界面
            for (ULL faceId : boundaryFaceIndexes)  // 直接用面上的phi值
            {
                matrix.addB(faces[faceId].getOwnerIndex(), -phiFaceField[faceId] * mf[faceId]);
            }
        }
        else if (type == DivType::MINMOD)
        {

        }


    }
}





#endif // DIV_H_