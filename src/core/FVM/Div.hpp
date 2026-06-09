#ifndef DIV_H_
#define DIV_H_

#include "Field.hpp"
#include "SparseMatrix.hpp"
#include "DivType.h"
#include "interpolation.hpp"
#include <cassert>


namespace fvm {

    namespace detail {

        /**
         * @brief 防止分母为 0 的安全除法
         */
        inline Scalar safeRatio(Scalar numerator, Scalar denominator)
        {
            constexpr Scalar EPSILON = 1e-14;

            if (std::abs(denominator) < EPSILON)
            {
                return 0.0;
            }

            return numerator / denominator;
        }


        /**
         * @brief MINMOD 限制器
         *
         * psi(r) = max(0, min(1, r))
         */
        inline Scalar minmodLimiter(Scalar r)
        {
            return std::max(0.0, std::min(1.0, r));
        }


        /**
         * @brief MUSCL 限制器
         *
         * psi(r) = max(0, min(2r, 0.5(1 + r), 2))
         */
        inline Scalar musclLimiter(Scalar r)
        {
            return std::max(
                0.0,
                std::min(
                    std::min(2.0 * r, 0.5 * (1.0 + r)),
                    2.0
                )
            );
        }


        /**
         * @brief 根据格式类型选择限制器
         */
        inline Scalar limiter(Scalar r, DivType type)
        {
            if (type == DivType::MINMOD)
            {
                return minmodLimiter(r);
            }

            if (type == DivType::MUSCL)
            {
                return musclLimiter(r);
            }

            return 0.0;
        }
    }

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
        else if (type == DivType::SUD)  // 二阶迎风：一阶迎风 + 二阶延迟修正
        {
            /*
             * 说明：
             *
             * 1. 当 Tp == Scalar 时：
             *    使用带 MINMOD 限制器的 SUD。
             *
             * 2. 当 Tp != Scalar 时：
             *    保持你原来的 SUD 写法。
             *
             * 原因：
             *    MINMOD/MUSCL 的 r = deltaUp / deltaDown 是标量比值。
             *    对 Vector/Tensor 场，不能直接定义这个标量 r。
             */
            if constexpr (std::is_same_v<Tp, Scalar>)
            {
                // ===============================
                // Scalar 场：带限制器的 SUD
                // ===============================

                // 先遍历内部面
                for (ULL faceId : internalFaceIndexes)
                {
                    const Face& face = faces[faceId];

                    const ULL owner = face.getOwnerIndex();
                    const ULL neighbor = face.getNeighborIndex();

                    const Scalar m_f = mf[faceId];

                    const Vector<Scalar> CD =
                        cells[neighbor].getCenter() - cells[owner].getCenter();

                    if (m_f >= 0.0)
                    {
                        /*
                         * 流向：
                         *
                         *     U -> owner(C) -> f -> neighbor(D)
                         *
                         * FUD:
                         *
                         *     phi_f^FUD = phi_C
                         *
                         * 你原来的 SUD 虚拟上上游构造：
                         *
                         *     phi_U = phi_D - 2 * grad_C · CD
                         *
                         * 限制器格式：
                         *
                         *     phi_f = phi_C + 0.5 * psi(r) * (phi_D - phi_C)
                         *
                         *     r = (phi_C - phi_U) / (phi_D - phi_C)
                         */

                         // 一阶迎风隐式部分
                        matrix(owner, owner) += m_f;
                        matrix(neighbor, owner) -= m_f;

                        const Scalar phiC = phiCellField_0[owner];
                        const Scalar phiD = phiCellField_0[neighbor];

                        const Scalar gradTerm =
                            cellGradient[owner] & CD;

                        const Scalar phiU =
                            phiD - 2.0 * gradTerm;

                        const Scalar deltaDown =
                            phiD - phiC;

                        const Scalar deltaUp =
                            phiC - phiU;

                        const Scalar r =
                            detail::safeRatio(deltaUp, deltaDown);

                        /*
                         * 先用 MINMOD，稳定优先。
                         * 如果之后觉得太耗散，再改成 detail::musclLimiter(r)。
                         */
                        const Scalar psi =
                            detail::minmodLimiter(r);

                        // const Scalar psi = detail::musclLimiter(r);

                        /*
                         * 受限制的二阶修正：
                         *
                         *     corr = phi_f - phi_f^FUD
                         *          = 0.5 * psi(r) * (phiD - phiC)
                         */
                        const Scalar corr =
                            0.5 * psi * deltaDown;

                        // 延迟修正源项
                        matrix.addB(owner, -m_f * corr);
                        matrix.addB(neighbor, m_f * corr);
                    }
                    else
                    {
                        /*
                         * 流向：
                         *
                         *     U -> neighbor(C) -> f -> owner(D)
                         *
                         * 此时：
                         *
                         *     C = neighbor
                         *     D = owner
                         *
                         * 注意：
                         *
                         *     CD = x_neighbor - x_owner
                         *
                         * 所以：
                         *
                         *     x_D - x_C = x_owner - x_neighbor = -CD
                         *
                         * 虚拟上上游构造：
                         *
                         *     phi_U = phi_D - 2 * grad_C · (x_D - x_C)
                         *           = phi_owner + 2 * grad_neighbor · CD
                         */

                         // 一阶迎风隐式部分
                        matrix(owner, neighbor) += m_f;
                        matrix(neighbor, neighbor) -= m_f;

                        const Scalar phiC = phiCellField_0[neighbor];
                        const Scalar phiD = phiCellField_0[owner];

                        const Scalar gradTerm =
                            cellGradient[neighbor] & CD;

                        const Scalar phiU =
                            phiD + 2.0 * gradTerm;

                        const Scalar deltaDown =
                            phiD - phiC;

                        const Scalar deltaUp =
                            phiC - phiU;

                        const Scalar r =
                            detail::safeRatio(deltaUp, deltaDown);

                        /*
                         * 先用 MINMOD。
                         */
                        const Scalar psi =
                            detail::minmodLimiter(r);

                        // const Scalar psi = detail::musclLimiter(r);

                        const Scalar corr =
                            0.5 * psi * deltaDown;

                        // 延迟修正源项
                        matrix.addB(owner, -m_f * corr);
                        matrix.addB(neighbor, m_f * corr);
                    }
                }

                // 再遍历边界面
                for (ULL faceId : boundaryFaceIndexes)
                {
                    const Face& face = faces[faceId];

                    const ULL owner = face.getOwnerIndex();
                    const Scalar m_f = mf[faceId];

                    if (m_f >= 0.0)
                    {
                        /*
                         * 边界出流：
                         *
                         *     内部 -> 边界
                         *
                         * 没有明确的下游单元 D 和上上游点 U，
                         * 所以这里保持退化为 FUD，更稳定。
                         */
                        matrix(owner, owner) += m_f;
                    }
                    else
                    {
                        /*
                         * 边界流入：
                         *
                         *     边界 -> 内部
                         *
                         * 直接使用边界面值。
                         */
                        matrix.addB(owner, -m_f * phiFaceField[faceId]);
                    }
                }
            }
            else
            {
                // ===============================
                // Vector / Tensor 场：保持你原来的 SUD 写法
                // ===============================

                // 先遍历内部面
                for (ULL faceId : internalFaceIndexes)
                {
                    const Face& face = faces[faceId];

                    const ULL owner = face.getOwnerIndex();
                    const LL neighbor = face.getNeighborIndex();

                    const Scalar m_f = mf[faceId];

                    const Vector<Scalar> CD =
                        cells[neighbor].getCenter() - cells[owner].getCenter();

                    const GradType ownerGrad = cellGradient[owner];
                    const GradType neighborGrad = cellGradient[neighbor];

                    if (m_f >= 0.0)
                    {
                        // 一阶部分
                        matrix(owner, owner) += m_f;
                        matrix(neighbor, owner) -= m_f;

                        // 二阶修正部分，保持你原来的泛型写法
                        matrix.addB(
                            owner,
                            -m_f * 1.5 * phiCellField_0[owner]
                            + m_f * 0.5 * phiCellField_0[neighbor]
                            - m_f * (ownerGrad & CD)
                            + m_f * phiCellField_0[owner]
                        );

                        matrix.addB(
                            neighbor,
                            -m_f * 0.5 * phiCellField_0[neighbor]
                            + m_f * 1.5 * phiCellField_0[owner]
                            + m_f * (ownerGrad & CD)
                            - m_f * phiCellField_0[owner]
                        );
                    }
                    else
                    {
                        // 一阶部分
                        matrix(owner, neighbor) += m_f;
                        matrix(neighbor, neighbor) -= m_f;

                        // 二阶修正部分，保持你原来的泛型写法
                        matrix.addB(
                            owner,
                            +m_f * 0.5 * phiCellField_0[owner]
                            - m_f * 1.5 * phiCellField_0[neighbor]
                            + m_f * (neighborGrad & CD)
                            + m_f * phiCellField_0[neighbor]
                        );

                        matrix.addB(
                            neighbor,
                            +m_f * 1.5 * phiCellField_0[neighbor]
                            - m_f * 0.5 * phiCellField_0[owner]
                            - m_f * (neighborGrad & CD)
                            - m_f * phiCellField_0[neighbor]
                        );
                    }
                }

                // 再遍历边界面
                for (ULL faceId : boundaryFaceIndexes)
                {
                    const Face& face = faces[faceId];

                    const ULL owner = face.getOwnerIndex();
                    const Scalar m_f = mf[faceId];

                    if (m_f >= 0.0)
                    {
                        // 一阶部分
                        matrix(owner, owner) += m_f;
                    }
                    else
                    {
                        matrix.addB(owner, -m_f * phiFaceField[faceId]);
                    }
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
        else if (type == DivType::MINMOD || type == DivType::MUSCL)
        {
            /*
             * MINMOD / MUSCL 目前只建议用于标量输运方程。
             *
             * 对 Vector / Tensor 场，r = deltaUp / deltaDown 不再是标量，
             * 需要分量限制器或者更复杂的多维限制器。
             */
            if constexpr (!std::is_same_v<Tp, Scalar>)
            {
                throw std::invalid_argument(
                    "Div() Error: MINMOD/MUSCL schemes are only implemented for Scalar fields."
                );
            }
            else
            {
                // ===============================
                // 1. 内部面
                // ===============================
                for (ULL faceId : internalFaceIndexes)
                {
                    const Face& face = faces[faceId];

                    const ULL owner = face.getOwnerIndex();
                    const ULL neighbor = face.getNeighborIndex();

                    const Scalar m_f = mf[faceId];

                    const Vector<Scalar> CD =
                        cells[neighbor].getCenter() - cells[owner].getCenter();

                    if (m_f >= 0.0)
                    {
                        /*
                         * 流向：
                         *
                         *     U ---- owner(C) ---- f ---- neighbor(D)
                         *
                         * FUD:
                         *
                         *     phi_f^FUD = phi_owner
                         *
                         * 虚拟上上游点：
                         *
                         *     phi_U ≈ phi_D - 2 * grad_C · CD
                         *
                         * TVD/MUSCL:
                         *
                         *     phi_f = phi_C + 0.5 * psi(r) * (phi_D - phi_C)
                         *
                         * 其中：
                         *
                         *     r = (phi_C - phi_U) / (phi_D - phi_C)
                         */

                         // 一阶迎风隐式部分
                        matrix(owner, owner) += m_f;
                        matrix(neighbor, owner) -= m_f;

                        const Scalar phiC = phiCellField_0[owner];
                        const Scalar phiD = phiCellField_0[neighbor];

                        const Scalar phiU =
                            phiD - 2.0 * (cellGradient[owner] & CD);

                        const Scalar deltaDown = phiD - phiC;
                        const Scalar deltaUp = phiC - phiU;

                        const Scalar r =
                            detail::safeRatio(deltaUp, deltaDown);

                        const Scalar psi =
                            detail::limiter(r, type);

                        /*
                         * corr = phi_f - phi_FUD
                         *      = 0.5 * psi(r) * (phi_D - phi_C)
                         */
                        const Scalar corr =
                            0.5 * psi * deltaDown;

                        // 延迟修正源项
                        matrix.addB(owner, -m_f * corr);
                        matrix.addB(neighbor, m_f * corr);
                    }
                    else
                    {
                        /*
                         * 流向：
                         *
                         *     U ---- neighbor(C) ---- f ---- owner(D)
                         *
                         * 注意：
                         *
                         *     CD = x_neighbor - x_owner
                         *
                         * 因此：
                         *
                         *     x_D - x_C = x_owner - x_neighbor = -CD
                         *
                         * 虚拟上上游点：
                         *
                         *     phi_U ≈ phi_D - 2 * grad_C · (x_D - x_C)
                         *           = phi_owner + 2 * grad_neighbor · CD
                         */

                         // 一阶迎风隐式部分
                        matrix(owner, neighbor) += m_f;
                        matrix(neighbor, neighbor) -= m_f;

                        const Scalar phiC = phiCellField_0[neighbor];
                        const Scalar phiD = phiCellField_0[owner];

                        const Scalar phiU =
                            phiD + 2.0 * (cellGradient[neighbor] & CD);

                        const Scalar deltaDown = phiD - phiC;
                        const Scalar deltaUp = phiC - phiU;

                        const Scalar r =
                            detail::safeRatio(deltaUp, deltaDown);

                        const Scalar psi =
                            detail::limiter(r, type);

                        /*
                         * corr = phi_f - phi_FUD
                         *      = 0.5 * psi(r) * (phi_D - phi_C)
                         */
                        const Scalar corr =
                            0.5 * psi * deltaDown;

                        // 延迟修正源项
                        matrix.addB(owner, -m_f * corr);
                        matrix.addB(neighbor, m_f * corr);
                    }
                }

                // ===============================
                // 2. 边界面
                // ===============================
                for (ULL faceId : boundaryFaceIndexes)
                {
                    const Face& face = faces[faceId];

                    const ULL owner = face.getOwnerIndex();
                    const Scalar m_f = mf[faceId];

                    if (m_f >= 0.0)
                    {
                        /*
                         * 边界出流：
                         *
                         *     内部 -> 边界
                         *
                         * 对 MINMOD / MUSCL，这里没有明确的下游单元 D，
                         * 所以最稳妥的处理是退化为 FUD。
                         *
                         * 这也是很多 TVD 格式在边界附近的常见处理。
                         */
                        matrix(owner, owner) += m_f;
                    }
                    else
                    {
                        /*
                         * 边界流入：
                         *
                         *     边界 -> 内部
                         *
                         * 直接使用边界面值。
                         */
                        matrix.addB(owner, -m_f * phiFaceField[faceId]);
                    }
                }
            }
        }
    }
}





#endif // DIV_H_