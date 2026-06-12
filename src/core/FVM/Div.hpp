#ifndef DIV_H_
#define DIV_H_

#include "Field.hpp"
#include "SparseMatrix.hpp"
#include "DivType.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <type_traits>


namespace fvm {
    using ULL = unsigned long long;
    namespace detail {

        /**
         * @brief 防止分母为 0 的安全除法。
         */
        inline Scalar safeRatio(const Scalar numerator, const Scalar denominator)
        {
            constexpr Scalar EPSILON = 1.0e-14;

            if (std::abs(denominator) < EPSILON)
            {
                return 0.0;
            }

            return numerator / denominator;
        }

        /**
         * @brief 将标量限制在 [lower, upper]。
         */
        inline Scalar clampScalar(
            const Scalar value,
            const Scalar lower,
            const Scalar upper
        )
        {
            return std::max(lower, std::min(value, upper));
        }

        /**
         * @brief MINMOD 限制器。
         *
         * psi(r) = max(0, min(1, r))
         */
        inline Scalar minmodLimiter(const Scalar r)
        {
            return std::max(0.0, std::min(1.0, r));
        }

        /**
         * @brief MUSCL / van Leer 族常用限制器形式。
         *
         * psi(r) = max(0, min(2r, 0.5(1 + r), 2))
         */
        inline Scalar musclLimiter(const Scalar r)
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
         * @brief 根据对流格式选择 TVD 限制器。
         */
        inline Scalar limiter(const Scalar r, const DivType type)
        {
            switch (type)
            {
            case DivType::MINMOD:
                return minmodLimiter(r);

            case DivType::MUSCL:
                return musclLimiter(r);

            default:
                return 1.0;
            }
        }

        /**
         * @brief 判断当前 DivType 是否为高阶格式。
         */
        inline bool isHighOrderScheme(const DivType type)
        {
            return (
                type == DivType::CD ||
                type == DivType::SUD ||
                type == DivType::MINMOD ||
                type == DivType::MUSCL
                );
        }

        /**
         * @brief 标量局部最小值。
         */
        inline Scalar componentMin(const Scalar lhs, const Scalar rhs)
        {
            return std::min(lhs, rhs);
        }

        /**
         * @brief 标量局部最大值。
         */
        inline Scalar componentMax(const Scalar lhs, const Scalar rhs)
        {
            return std::max(lhs, rhs);
        }

        /**
         * @brief 矢量逐分量局部最小值。
         */
        inline Vector<Scalar> componentMin(
            const Vector<Scalar>& lhs,
            const Vector<Scalar>& rhs
        )
        {
            return Vector<Scalar>(
                std::min(lhs.x(), rhs.x()),
                std::min(lhs.y(), rhs.y()),
                std::min(lhs.z(), rhs.z())
            );
        }

        /**
         * @brief 矢量逐分量局部最大值。
         */
        inline Vector<Scalar> componentMax(
            const Vector<Scalar>& lhs,
            const Vector<Scalar>& rhs
        )
        {
            return Vector<Scalar>(
                std::max(lhs.x(), rhs.x()),
                std::max(lhs.y(), rhs.y()),
                std::max(lhs.z(), rhs.z())
            );
        }

        /**
         * @brief 标量有界限制。
         */
        inline Scalar clampValue(
            const Scalar value,
            const Scalar lower,
            const Scalar upper
        )
        {
            return clampScalar(value, lower, upper);
        }

        /**
         * @brief 矢量逐分量有界限制。
         */
        inline Vector<Scalar> clampValue(
            const Vector<Scalar>& value,
            const Vector<Scalar>& lower,
            const Vector<Scalar>& upper
        )
        {
            return Vector<Scalar>(
                clampScalar(value.x(), lower.x(), upper.x()),
                clampScalar(value.y(), lower.y(), upper.y()),
                clampScalar(value.z(), lower.z(), upper.z())
            );
        }

        /**
         * @brief 标量 TVD 限制。
         *
         * rawCorrection 是从上游单元中心到面心的高阶修正：
         *
         *     phi_f^HO = phi_U + rawCorrection
         *
         * downwindJump 是：
         *
         *     phi_D - phi_U
         *
         * 对于均匀一维网格，rawCorrection ≈ 0.5 * (phi_D - phi_U)。
         * 因此这里用 0.5 * downwindJump 作为归一化尺度构造 r。
         */
        inline Scalar applyTvdLimiter(
            const Scalar rawCorrection,
            const Scalar downwindJump,
            const DivType type
        )
        {
            if (type != DivType::MINMOD && type != DivType::MUSCL)
            {
                return rawCorrection;
            }

            const Scalar referenceCorrection = 0.5 * downwindJump;
            const Scalar r = safeRatio(rawCorrection, referenceCorrection);
            const Scalar psi = limiter(r, type);

            return psi * rawCorrection;
        }

        /**
         * @brief 矢量逐分量 TVD 限制。
         */
        inline Vector<Scalar> applyTvdLimiter(
            const Vector<Scalar>& rawCorrection,
            const Vector<Scalar>& downwindJump,
            const DivType type
        )
        {
            if (type != DivType::MINMOD && type != DivType::MUSCL)
            {
                return rawCorrection;
            }

            return Vector<Scalar>(
                applyTvdLimiter(rawCorrection.x(), downwindJump.x(), type),
                applyTvdLimiter(rawCorrection.y(), downwindJump.y(), type),
                applyTvdLimiter(rawCorrection.z(), downwindJump.z(), type)
            );
        }

        /**
         * @brief 沿 owner-neighbour 中心连线计算面心线性插值权重。
         *
         * 非结构网格上不建议只用法向投影距离，因为：
         *
         *     C_N - C_P 不一定平行于 S_f
         *
         * 这里采用：
         *
         *     lambda = ((C_f - C_P) · (C_N - C_P)) / |C_N - C_P|^2
         *
         * 然后限制到 [0, 1]。
         */
        inline Scalar centreLineWeight(
            const Vector<Scalar>& ownerCenter,
            const Vector<Scalar>& neighborCenter,
            const Vector<Scalar>& faceCenter
        )
        {
            const Vector<Scalar> ownerToNeighbor = neighborCenter - ownerCenter;
            const Vector<Scalar> ownerToFace = faceCenter - ownerCenter;

            const Scalar denominator = ownerToNeighbor & ownerToNeighbor;

            if (std::abs(denominator) < 1.0e-30)
            {
                return 0.5;
            }

            return clampScalar(
                (ownerToFace & ownerToNeighbor) / denominator,
                0.0,
                1.0
            );
        }

        /**
         * @brief 线性插值。
         *
         * lambda = 0 -> ownerValue
         * lambda = 1 -> neighborValue
         */
        template<typename Tp>
        inline Tp linearInterpolate(
            const Tp& ownerValue,
            const Tp& neighborValue,
            const Scalar lambda
        )
        {
            return (1.0 - lambda) * ownerValue + lambda * neighborValue;
        }

        /**
         * @brief 计算上游梯度到面心的重构修正。
         *
         * Scalar 场：
         *
         *     GradType = Vector<Scalar>
         *     gradPhi & d -> Scalar
         *
         * Vector 场：
         *
         *     GradType = Tensor<Scalar>
         *     gradU & d -> Vector<Scalar>
         */
        template<typename GradType>
        inline auto gradientCorrection(
            const GradType& grad,
            const Vector<Scalar>& displacement
        ) -> decltype(grad& displacement)
        {
            return grad & displacement;
        }

        /**
         * @brief 计算一个内部面的 bounded deferred correction。
         *
         * 所有高阶格式都写成：
         *
         *     phi_f = phi_upwind + correction
         *
         * 其中 correction 显式加入右端项。
         *
         * 为避免非结构网格上的高阶外推振荡，最终面值强制限制在
         *
         *     [min(phi_owner, phi_neighbor), max(phi_owner, phi_neighbor)]
         *
         * 对 Vector<Scalar> 是逐分量限制。
         */
        template<typename Tp, typename GradType>
        inline Tp boundedInternalCorrection(
            const DivType type,
            const Face& face,
            const Cell& ownerCell,
            const Cell& neighborCell,
            const Tp& ownerValue,
            const Tp& neighborValue,
            const GradType& ownerGradient,
            const GradType& neighborGradient,
            const Scalar massFlux
        )
        {
            const bool ownerIsUpwind = (massFlux >= 0.0);

            const Tp& upwindValue =
                ownerIsUpwind ? ownerValue : neighborValue;

            const Tp& downwindValue =
                ownerIsUpwind ? neighborValue : ownerValue;

            const Vector<Scalar>& upwindCenter =
                ownerIsUpwind ? ownerCell.getCenter() : neighborCell.getCenter();

            const GradType& upwindGradient =
                ownerIsUpwind ? ownerGradient : neighborGradient;

            Tp candidateFaceValue{};

            if (type == DivType::CD)
            {
                /*
                 * bounded CD:
                 *
                 * 先计算几何中心线插值值，然后再限制。
                 * 裸中心差分在强对流或间断问题中不具备单调性，
                 * 所以这里不使用纯隐式 CD 矩阵，而是采用 FUD + bounded correction。
                 */
                const Scalar lambda = centreLineWeight(
                    ownerCell.getCenter(),
                    neighborCell.getCenter(),
                    face.getCenter()
                );

                candidateFaceValue = linearInterpolate(
                    ownerValue,
                    neighborValue,
                    lambda
                );
            }
            else
            {
                /*
                 * SUD / MINMOD / MUSCL:
                 *
                 * 非结构网格通用写法：
                 *
                 *     phi_f^HO = phi_U + grad(phi)_U · (C_f - C_U)
                 *
                 * 对 MINMOD / MUSCL，再对 correction 做 TVD 限制。
                 */
                const Vector<Scalar> displacement =
                    face.getCenter() - upwindCenter;

                Tp rawCorrection =
                    gradientCorrection(upwindGradient, displacement);

                rawCorrection = applyTvdLimiter(
                    rawCorrection,
                    downwindValue - upwindValue,
                    type
                );

                candidateFaceValue = upwindValue + rawCorrection;
            }

            /*
             * 局部有界限制。
             *
             * 注意：
             * 这一步是抑制非结构网格振荡的关键。
             */
            const Tp lower = componentMin(ownerValue, neighborValue);
            const Tp upper = componentMax(ownerValue, neighborValue);

            const Tp boundedFaceValue =
                clampValue(candidateFaceValue, lower, upper);

            return boundedFaceValue - upwindValue;
        }

        /**
         * @brief 添加内部面的 FUD 隐式矩阵贡献。
         */
        template<typename Tp>
        inline void addImplicitFudInternalFace(
            SparseMatrix<Tp>& matrix,
            const ULL owner,
            const ULL neighbor,
            const Scalar massFlux
        )
        {
            if (massFlux >= 0.0)
            {
                matrix(owner, owner) += massFlux;
                matrix(neighbor, owner) -= massFlux;
            }
            else
            {
                matrix(owner, neighbor) += massFlux;
                matrix(neighbor, neighbor) -= massFlux;
            }
        }

        /**
         * @brief 添加内部面的 deferred correction 右端项。
         *
         * owner 方程中，当前面的通量贡献符号为 +mf。
         * neighbor 方程中，当前面的通量贡献符号为 -mf。
         *
         * 若：
         *
         *     phi_f = phi_upwind + correction
         *
         * 则：
         *
         *     owner RHS   += -mf * correction
         *     neighbor RHS +=  mf * correction
         */
        template<typename Tp>
        inline void addDeferredCorrectionInternalFace(
            SparseMatrix<Tp>& matrix,
            const ULL owner,
            const ULL neighbor,
            const Scalar massFlux,
            const Tp& correction
        )
        {
            matrix.addB(owner, -massFlux * correction);
            matrix.addB(neighbor, massFlux * correction);
        }

        /**
         * @brief 边界面对流项处理。
         *
         * 对边界面，统一采用：
         *
         * 1. 出流 mf >= 0：
         *
         *        phi_f = phi_owner
         *
         * 2. 入流 mf < 0：
         *
         *        phi_f = phi_boundaryFace
         *
         * 高阶格式在边界附近没有完整 owner-neighbor 模板，
         * 强行外推容易产生边界振荡，因此这里稳定优先。
         */
        template<typename Tp>
        inline void addBoundaryFaceContribution(
            SparseMatrix<Tp>& matrix,
            const ULL owner,
            const Scalar massFlux,
            const Tp& boundaryFaceValue
        )
        {
            if (massFlux >= 0.0)
            {
                matrix(owner, owner) += massFlux;
            }
            else
            {
                matrix.addB(owner, -massFlux * boundaryFaceValue);
            }
        }

    } // namespace detail


    /**
     * @brief 对流项离散函数。
     *
     * 支持：
     *
     *     Scalar
     *     Vector<Scalar>
     *
     * 离散形式：
     *
     *     div(rho U phi)
     *
     * 内部面通量：
     *
     *     mf = rho_f * (U_f · S_f)
     *
     * FUD：
     *
     *     完全隐式一阶迎风。
     *
     * CD / SUD / MINMOD / MUSCL：
     *
     *     FUD 隐式矩阵 + bounded deferred correction。
     */
    template<typename Tp>
    void Div(
        SparseMatrix<Tp>& matrix,
        const FaceField<Scalar>& rho,
        Field<Tp>& phi,
        FaceField<Vector<Scalar>>& U,
        DivType type
    );


    template<typename Tp>
    void Div(
        SparseMatrix<Tp>& matrix,
        const FaceField<Scalar>& rho,
        Field<Tp>& phi,
        FaceField<Vector<Scalar>>& U,
        DivType type
    )
    {
        static_assert(
            std::is_same_v<Tp, Scalar> ||
            std::is_same_v<Tp, Vector<Scalar>>,
            "fvm::Div() only supports Scalar and Vector<Scalar> fields."
            );

        using ULL = unsigned long long;
        using LL = long long;
        using GradType = decltype(Tp()* Vector<Scalar>());

        Mesh* mesh = matrix.getMesh();

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

        if (!U.isValid())
        {
            std::cerr << "Div() Error: U is not valid." << std::endl;
            throw std::invalid_argument("U is not valid.");
        }

        if (mesh == nullptr)
        {
            std::cerr << "Div() Error: mesh is nullptr." << std::endl;
            throw std::invalid_argument("mesh is nullptr.");
        }

        if (mesh != rho.getMesh() ||
            mesh != phi.getMesh() ||
            mesh != U.getMesh())
        {
            std::cerr
                << "Div() Error: The mesh of matrix, rho, phi and U is not the same."
                << std::endl;

            throw std::invalid_argument(
                "The mesh of matrix, rho, phi and U is not the same."
            );
        }

        if (matrix.fieldPtr_ == nullptr)
        {
            matrix.fieldPtr_ = &phi;
        }
        else if (matrix.fieldPtr_ != &phi)
        {
            std::cerr
                << "Div() Error: The field of the matrix is not the same as phi."
                << std::endl;

            throw std::invalid_argument(
                "The field of the matrix is not the same as phi."
            );
        }

        if (!detail::isHighOrderScheme(type) && type != DivType::FUD)
        {
            std::cerr << "Div() Error: unknown DivType." << std::endl;
            throw std::invalid_argument("unknown DivType.");
        }

        const std::vector<Face>& faces = mesh->getFaces();
        const std::vector<Cell>& cells = mesh->getCells();

        const std::vector<ULL>& internalFaceIndexes =
            mesh->getInternalFaceIndexes();

        const std::vector<ULL>& boundaryFaceIndexes =
            mesh->getBoundaryFaceIndexes();

        const FaceField<Tp>& phiFaceField =
            phi.getFaceField();

        const CellField<Tp>& phiCellField =
            phi.getCellField();

        /*
         * 高阶格式需要梯度。
         *
         * 你的 Field::cellToFace() 会更新 cellGradientField_。
         * 因此调用高阶格式前，外部应先调用：
         *
         *     phi.cellToFace();
         */
        const CellField<GradType>& cellGradient =
            phi.getCellGradientField();

        /*
         * 1. 预计算所有面的质量流量：
         *
         *     mf = rho_f * U_f · S_f
         *
         * 这里 face.getNormal() 已经在 Mesh::calculateMeshInfo()
         * 中修正为从 owner 指向 neighbor / 外部。
         */
        std::vector<Scalar> massFlux(mesh->getFaceNumber(), 0.0);

        for (ULL faceId = 0; faceId < mesh->getFaceNumber(); ++faceId)
        {
            const Face& face = faces[faceId];

            massFlux[faceId] =
                rho[faceId] *
                (U[faceId] & (face.getNormal() * face.getArea()));
        }

        /*
         * 2. 内部面。
         */
        for (const ULL faceId : internalFaceIndexes)
        {
            const Face& face = faces[faceId];

            const ULL owner = face.getOwnerIndex();
            const LL neighborSigned = face.getNeighborIndex();

            if (neighborSigned < 0)
            {
                std::cerr
                    << "Div() Error: internal face has invalid neighbor index."
                    << std::endl;

                throw std::runtime_error(
                    "internal face has invalid neighbor index."
                );
            }

            const ULL neighbor = static_cast<ULL>(neighborSigned);
            const Scalar mf = massFlux[faceId];

            /*
             * 所有格式都先添加 FUD 隐式部分。
             *
             * 这样可以保证矩阵主体是上风型，避免纯高阶格式破坏矩阵单调性。
             */
            detail::addImplicitFudInternalFace(
                matrix,
                owner,
                neighbor,
                mf
            );

            if (type == DivType::FUD)
            {
                continue;
            }

            /*
             * 高阶格式：
             *
             *     bounded correction = phi_f^bounded - phi_upwind
             *
             * 然后将 correction 作为显式源项加入右端项。
             */
            const Tp correction =
                detail::boundedInternalCorrection(
                    type,
                    face,
                    cells[owner],
                    cells[neighbor],
                    phiCellField[owner],
                    phiCellField[neighbor],
                    cellGradient[owner],
                    cellGradient[neighbor],
                    mf
                );

            detail::addDeferredCorrectionInternalFace(
                matrix,
                owner,
                neighbor,
                mf,
                correction
            );
        }

        /*
         * 3. 边界面。
         *
         * 边界统一采用稳定处理：
         *
         *     出流：owner 上风
         *     入流：边界面值
         */
        for (const ULL faceId : boundaryFaceIndexes)
        {
            const Face& face = faces[faceId];

            const ULL owner = face.getOwnerIndex();
            const Scalar mf = massFlux[faceId];

            detail::addBoundaryFaceContribution(
                matrix,
                owner,
                mf,
                phiFaceField[faceId]
            );
        }
    }

} // namespace fvm


#endif // DIV_H_