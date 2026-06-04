#ifndef SIMPLE1_H_
#define SIMPLE1_H_

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

#include "Div.hpp"
#include "DivType.h"
#include "Field.hpp"
#include "Laplacian.hpp"
#include "Mesh.h"
#include "Solver.hpp"
#include "SparseMatrix.hpp"
#include "Vector.hpp"
#include "Gradient.hpp"
#include "Divergence.hpp"
#include "interpolation.hpp"

namespace algorithm {
    namespace simple {
        class SIMPLE {
            using ULL = unsigned long long;
            using LL = long long;
            using Scalar = double;
        public:
            /**
             * @brief SIMPLE算法参数设置
             */
            struct Options {
                // SIMPLE 外迭代控制
                int maxOuterIterations = 500;
                Scalar convergenceTolerance = 1e-8;
                int logInterval = 20;

                // 第二类边界条件内迭代
                int innerMaxIterations = 5000;
                Scalar relativeTolerance = 1e-6;

                // 动量方程内部迭代
                int momentumMaxIterations = 5000;
                Scalar momentumTolerance = 1e-10;
                Solver<Vector<Scalar>>::Method momentumMethod =
                    Solver<Vector<Scalar>>::Method::Jacobi;
                Solver<Scalar>::Method pressureMethod =
                    Solver<Scalar>::Method::Jacobi;

                // 压力修正方程内部迭代
                int pressureMaxIterations = 10000;
                Scalar pressureTolerance = 1e-10;

                // SIMPLE 欠松弛
                Scalar alphaU = 0.7;
                Scalar alphaP = 0.3;

                // 压力参考
                ULL pressureReferenceCell = 0;
                Scalar pressureReferenceValue = 0.0;
                Scalar pressureReferencePenalty = 1e20;

                // 非正交修正次数。
                // 1 表示只做正交压力修正；大于 1 时会使用显式非正交修正。
                int nNonOrthogonalCorrectors = 1;

                // 对流格式
                fvm::DivType divScheme = fvm::DivType::FUD;

                // 是否使用 Rhie-Chow 修正面通量
                bool useRhieChow = true;


                // 压力修正边界 p'=0 的 patch 名称。
                // 比如压力出口可填 {"outlet"}。
                // 未列出的边界默认按 dp'/dn=0 处理。
                std::vector<std::string> fixedPressurePatches;
            };

            /**
             * @brief SIMPLE 算法残差收敛信息
             */
            struct Residual {
                int iteration = 0;
                Scalar continuity = std::numeric_limits<Scalar>::max();
                Scalar velocity = std::numeric_limits<Scalar>::max();
                Scalar pressure = std::numeric_limits<Scalar>::max();
                // Scalar pressureEquation = std::numeric_limits<Scalar>::max();

                bool converged(Scalar tolerance) const {
                    return continuity < tolerance &&
                        velocity < tolerance &&
                        pressure < tolerance;
                }
            };

            /**
             * @brief 构造函数，仅支持有参构造
             * @param U 初始速度场
             * @param p 初始压力场
             * @param rho 密度
             * @param mu 粘度
             */
            SIMPLE(Field<Vector<Scalar>>& U,
                   Field<Scalar>& p,
                   FaceField<Scalar>& rho,
                   FaceField<Scalar>& mu);

            /**
             * @brief 构造函数，仅支持有参构造
             * @param U 初始速度场
             * @param p 初始压力场
             * @param rho 密度
             * @param mu 粘度
             * @param options SIMPLE 算法参数
             */
            SIMPLE(Field<Vector<Scalar>>& U,
                   Field<Scalar>& p,
                   FaceField<Scalar>& rho,
                   FaceField<Scalar>& mu,
                   const Options& options);

            SIMPLE() = delete;
            SIMPLE(const SIMPLE&) = delete;
            SIMPLE(SIMPLE&&) = delete;
            SIMPLE& operator=(const SIMPLE&) = delete;
            SIMPLE& operator=(SIMPLE&&) = delete;

        public:
            /**
             * @brief 运行 SIMPLE 算法
             */
            void solve();

        private:
            /**
             * @brief 检查输入场是否合法
             */
            void isValidationInput();

            /**
             * @brief 检查算法参数是否合法
             */
            void isValidationOptions();

            /**
             * @brief 初始化 rAU_, 设置初始场值，并设置边界条件
             */
            void initializeRAU();

            /**
             * @brief 线性插值系数
             * @param faceId 面ID
             * @return
             */
            Scalar interpolationAlphaInternal(ULL faceId) const;

            /**
             * @brief grad(p)_f^face & Sf 计算辅助函数
             * @param internalFaceId 面ID
             * @param pOwner 主单元压力
             * @param pNeighbor 邻单元压力
             * @param gradPF 面压力梯度
             * @return
             */
            Scalar faceNormalPressureGradientDotSf(ULL internalFaceId,
                                                   Scalar pOwner,
                                                   Scalar pNeighbor,
                                                   const Vector<Scalar>& gradPF) const;
            /**
             * @brief 通过Rhie-Chow插值计算面速度Uf_
             */
            void computeRhieChowFaceVelocity();

            /**
             * @brief 计算动量方程矩阵
             */
            void solverMomentumEqn();

            /**
             * @brief 计算压力修正方程矩阵
             */
            void solverPressureCorrEqn();

            /**
             * @brief 计算 rAU_ = V / a_P
             */
            void computeRAU();

            /**
             * @brief 初始化压力修正场 pCorr_
             */
            void initializepCorr();

            /**
             * @brief 初始化固定压力面和边界条件
             */
            void initializeFixedPressureFaces();

            /**
             * @brief 修正计算压力修正值后修正Uf_
             */
            void correctFaceVelocity();

            /**
             * @brief 修正压力
             */
            void correctPressure();

            /**
             * @brief 修正速度
             */
            void correctVelocity();

            /**
             * @brief 计算质量守恒
             */
            Scalar computeContinuityResidual() const;

            /**
             * @brief 计算压力残差
             * @param current 本部场
             * @param old 上一步场
             * @return 残差
             */
            Scalar computeRelativeChange(const std::vector<Scalar>& current,
                                         const std::vector<Scalar>& old) const;

            /**
             * @brief 计算速度残差
             * @param current 本部场
             * @param old 上一步场
             * @return 残差
             */
            Scalar computeRelativeChange(const std::vector<Vector<Scalar>>& current,
                                         const std::vector<Vector<Scalar>>& old) const;

        private:
            static constexpr Scalar SMALL = 1e-30;

            // SIMPLE 算法参数
            Options options_;

            // 各种场
            Field<Vector<Scalar>>& U_;      // 速度场
            Field<Scalar>& p_;              // 压力场
            FaceField<Scalar>& rho_;        // 密度
            FaceField<Scalar>& mu_;         // 粘度

            Mesh* mesh_;                    // 网格

            Field<Scalar> rAU_;             // rAU = V / a_P   压力修正方程扩散系数
            FaceField<Vector<Scalar>> Uf_;  // 面场速度, 使用Rhie-Chow才启用
            Field<Scalar> pCorr_;           // 压力修正场

            // 方程组
            SparseMatrix<Vector<Scalar>> momentumEqn_;  // 动量方程矩阵
            SparseMatrix<Scalar> pressureCorrEqn_;      // 压力修正方程矩阵

            std::vector<bool> fixedPressureFace_;       // 固定压力面
        };  // class SIMPLE




    #pragma region 实现
        inline SIMPLE::SIMPLE(Field<Vector<Scalar>>& U,
                              Field<Scalar>& p,
                              FaceField<Scalar>& rho,
                              FaceField<Scalar>& mu)
            : options_(SIMPLE::Options())
            , U_(U)
            , p_(p)
            , rho_(rho)
            , mu_(mu)
            , mesh_(U.getMesh())
            , rAU_("rAU", mesh_)
            , Uf_("Uf", mesh_)
            , pCorr_("pCorr", mesh_)
            , momentumEqn_(mesh_)
            , pressureCorrEqn_(mesh_)
            , fixedPressureFace_(mesh_->getFaceNumber(), false)
        {
            // 检查参数是否合理
            isValidationInput();

            // 初始化 rAU_, pCorr_
            initializeRAU();
            initializepCorr();
        }

        inline SIMPLE::SIMPLE(Field<Vector<Scalar>>& U,
                              Field<Scalar>& p,
                              FaceField<Scalar>& rho,
                              FaceField<Scalar>& mu,
                              const Options& options)
            : options_(options)
            , U_(U)
            , p_(p)
            , rho_(rho)
            , mu_(mu)
            , mesh_(U.getMesh())
            , rAU_("rAU", mesh_)
            , Uf_("Uf", mesh_)
            , pCorr_("pCorr", mesh_)
            , momentumEqn_(mesh_)
            , pressureCorrEqn_(mesh_)
            , fixedPressureFace_(mesh_->getFaceNumber(), false)
        {
            // 检查参数是否合理
            isValidationInput();
            // 检查算法参数是否合理
            isValidationOptions();

            // 初始化 rAU_, pCorr_
            initializeRAU();
            initializepCorr();
        }
    #pragma region solver
        inline void SIMPLE::solve() {
            // 先插值 速度/压力 面场
            U_.cellToFace();
            p_.cellToFace();

            // 采用Rhie-Chow插值计算面场速度Uf_
            computeRhieChowFaceVelocity();

            Residual residual;
            // 外迭代
            for (int outerIter = 0; outerIter < options_.maxOuterIterations; ++outerIter) {
                residual.iteration = outerIter + 1;
                // 保留旧值
                std::vector<Scalar> oldP = p_.getCellField().getData();
                std::vector<Vector<Scalar>> oldU = U_.getCellField().getData();

                // 计算动量方程
                solverMomentumEqn();
                U_.cellToFace();
                // 计算Df
                computeRAU();

                computeRhieChowFaceVelocity();

                // 计算压力修正方程
                solverPressureCorrEqn();
                pCorr_.cellToFace();

                // 用 p' 修正 Uf*
                correctFaceVelocity();

                // 修正压力, 速度
                correctPressure();
                correctVelocity();

                // Rhie-Chow插值

                // 计算残差
                residual.continuity = computeContinuityResidual();
                residual.velocity = computeRelativeChange(
                    U_.getCellField().getData(),
                    oldU
                );
                residual.pressure = computeRelativeChange(
                    p_.getCellField().getData(),
                    oldP
                );

                // 打印输出
                if (options_.logInterval > 0 &&
                    ((outerIter + 1) % options_.logInterval == 0||
                     outerIter == 0 ||
                     residual.converged(options_.convergenceTolerance)))
                {
                    std::cout
                        << "SIMPLE outer iteration " << outerIter + 1
                        << ", continuity = " << residual.continuity
                        << ", U = " << residual.velocity
                        << ", p = " << residual.pressure
                        // << ", pEqn = " << residual.pressureEquation
                        << std::endl;
                }
                if (residual.converged(options_.convergenceTolerance)) {
                    break;
                }
            }
        }


        inline void SIMPLE::isValidationInput() {
            if (mesh_ == nullptr) {
                std::cerr << "SIMPLE::SIMPLE: mesh is nullptr" << std::endl;
                throw std::invalid_argument("mesh is nullptr");
            }
            if (mesh_ != rho_.getMesh() ||
                mesh_ != mu_.getMesh() ||
                mesh_ != p_.getMesh())
            {
                std::cerr << "SIMPLE::SIMPLE: mesh is not equal" << std::endl;
                throw std::invalid_argument("U_, rho_, mu_, p_ mesh must be equal");
            }

            if (!U_.isValid()) {
                std::cerr << "SIMPLE::SIMPLE: U_ is not valid" << std::endl;
                throw std::invalid_argument("U_ is not valid");
            }
            if (!p_.isValid()) {
                std::cerr << "SIMPLE::SIMPLE: p_ is not valid" << std::endl;
                throw std::invalid_argument("p_ is not valid");
            }
            if (!rho_.isValid()) {
                std::cerr << "SIMPLE::SIMPLE: rho_ is not valid" << std::endl;
                throw std::invalid_argument("rho_ is not valid");
            }
            if (!mu_.isValid()) {
                std::cerr << "SIMPLE::SIMPLE: mu_ is not valid" << std::endl;
                throw std::invalid_argument("mu_ is not valid");
            }
        }
        inline void SIMPLE::isValidationOptions() {
            if (options_.maxOuterIterations <= 0) {
                throw std::runtime_error(
                    "SIMPLE Error: maxOuterIterations must be positive."
                );
            }

            if (options_.momentumMaxIterations <= 0 ||
                options_.pressureMaxIterations <= 0) {
                throw std::runtime_error(
                    "SIMPLE Error: inner solver iteration numbers must be positive."
                );
            }

            if (options_.alphaU <= 0.0 || options_.alphaU > 1.0) {
                throw std::runtime_error(
                    "SIMPLE Error: alphaU must be in (0, 1]."
                );
            }

            if (options_.alphaP <= 0.0 || options_.alphaP > 1.0) {
                throw std::runtime_error(
                    "SIMPLE Error: alphaP must be in (0, 1]."
                );
            }

            if (options_.pressureReferenceCell >= mesh_->getCellNumber()) {
                throw std::runtime_error(
                    "SIMPLE Error: pressureReferenceCell is out of range."
                );
            }

            if (options_.nNonOrthogonalCorrectors <= 0) {
                throw std::runtime_error(
                    "SIMPLE Error: nNonOrthogonalCorrectors must be positive."
                );
            }
        }
        inline void SIMPLE::initializeRAU() {
            rAU_.setValue(Scalar());
            // 遍历边界，设置边界条件
            for (const auto& [name, patch] :
                 mesh_->getBoundaryPatches())
            {
                // 跳过二维EMPTY边界
                if (patch.getType() ==
                    BoundaryPatch::BoundaryType::EMPTY) {
                    continue;
                }
                rAU_.setBoundaryCondition(name, 0, 1, 0);
            }
            rAU_.cellToFace();
        }
        inline Scalar SIMPLE::interpolationAlphaInternal(ULL faceId) const
        {
            const std::vector<Face>& faces = mesh_->getFaces();
            const std::vector<Cell>& cells = mesh_->getCells();

            const Face& face = faces[faceId];

            const ULL owner = face.getOwnerIndex();
            const LL neighbor = face.getNeighborIndex();

            if (neighbor < 0) {
                std::cerr << "SIMPLE::SIMPLE: neighbor is negative" << std::endl;
                throw std::runtime_error("SIMPLE::SIMPLE: faceId must be internal faceId");
            }

            const Vector<Scalar>& faceCenter = face.getCenter();
            const Vector<Scalar>& ownerCenter = cells[owner].getCenter();
            const Vector<Scalar>& neighborCenter =
                cells[static_cast<ULL>(neighbor)].getCenter();

            const Vector<Scalar>& normal = face.getNormal();

            const Scalar ownerDistance =
                (faceCenter - ownerCenter) & normal;
            const Scalar neighborDistance =
                (neighborCenter - faceCenter) & normal;

            const Scalar denominator =
                ownerDistance + neighborDistance;

            Scalar alpha = ownerDistance / denominator;

            return alpha;
        }
        inline Scalar SIMPLE::faceNormalPressureGradientDotSf(ULL internalFaceId, Scalar pOwner, Scalar pNeighbor, const Vector<Scalar>& gradPF) const {
            const std::vector<Face>& faces = mesh_->getFaces();
            const std::vector<Cell>& cells = mesh_->getCells();

            const Face& face = faces[internalFaceId];
            const ULL owner = face.getOwnerIndex();
            const LL neighbor = face.getNeighborIndex();

            const Vector<Scalar> Sf = face.getNormal() * face.getArea();
            const Vector<Scalar> d_PN = cells[neighbor].getCenter() - cells[owner].getCenter();
            const Scalar distance = d_PN.magnitude();
            const Vector<Scalar> e = d_PN / distance;
            const Scalar denom = e & face.getNormal();
            const Scalar Efmag = face.getArea() / denom;

            const Vector<Scalar> Ef = Efmag * e;
            const Vector<Scalar> Tf = Sf - Ef;

            return (pNeighbor - pOwner) / distance * Efmag + (gradPF & Tf);
        }
        inline void SIMPLE::computeRhieChowFaceVelocity() {
            Uf_.setValue(Vector<Scalar>()); // 初始化为0, 并启用场(isValid)
            Uf_.getData() = U_.getFaceField().getData(); // 先初始化为 U_ 的面场值

            const std::vector<Face>& faces = mesh_->getFaces();
            const FaceField<Vector<Scalar>>& Uf = U_.getFaceField();
            const CellField<Scalar>& pCell = p_.getCellField();
            const CellField<Vector<Scalar>>& gradPCell = p_.getCellGradientField();
            const CellField<Scalar>& rAUp = rAU_.getCellField();
            const FaceField<Scalar>& rAUf = rAU_.getFaceField();

            // 通过Rhie-Chow插值计算内部面的Uf_场值
            for (ULL faceId : mesh_->getInternalFaceIndexes()) {
                const Face& face = faces[faceId];
                const ULL owner = face.getOwnerIndex();
                const LL neighbor = face.getNeighborIndex();

                const Vector<Scalar> Sf = face.getArea() * face.getNormal();
                const Scalar magSf2 = Sf.magnitudeSquared();

                // 插值权重：phi_f = (1-alpha) * phi_owner + alpha * phi_neighbor
                const Scalar alpha = interpolationAlphaInternal(faceId);

                // 插值得到面压力梯度
                Vector<Scalar> gradPf =
                    util::interpolate(gradPCell[owner],
                                      gradPCell[neighbor],
                                      interpolation::Scheme::LINEAR,
                                      alpha);

                // 面上的 D_f = V/a 插值
                const Scalar Df =
                    util::interpolate(rAUp[owner],
                                      rAUp[neighbor],
                                      interpolation::Scheme::LINEAR,
                                      alpha);

                // (grad p)_f^face · S_f 带非正交修正
                const Scalar gradPFaceDotSf =
                    faceNormalPressureGradientDotSf(faceId,
                                                    pCell[owner],
                                                    pCell[neighbor],
                                                    gradPf);

                // 体积流量修正 D_f * [gradPbar_f · S_f - (grad p)_f^face · S_f]
                const Scalar VolumeFluxCorr =
                    Df * ((gradPf & Sf) - gradPFaceDotSf);

                // 速度修正量只沿 S_f 方向修正，即只改法向速度。
                const Vector<Scalar> velocityCorrection =
                    (VolumeFluxCorr / magSf2) * Sf;

                // Uf = Uf + velocityCorrection
                const Vector<Scalar> URhieChow =
                    Uf[faceId] + velocityCorrection;
                Uf_[faceId] = URhieChow;
            }

            // 边界面
            for (ULL faceId : mesh_->getBoundaryFaceIndexes()) {
                Uf_[faceId] = Uf[faceId];
            }
        }
        inline void SIMPLE::solverMomentumEqn() {
            // 负压力梯度
            CellField<Vector<Scalar>> negGradPCell("-gradP", mesh_, Vector<Scalar>());
            std::vector<Vector<Scalar>>& negGradPCellData = negGradPCell.getData();
            for (std::size_t cellId = 0; const Vector<Scalar>& gradP : p_.getCellGradientField().getData()) {
                negGradPCellData[cellId] = -gradP;
                ++cellId;
            }
            // 内迭代(适应第二类边界条件)
            for (int i = 0; i < 1; ++i) {
                momentumEqn_.clear();
                // 对流项，扩散项，源项
                fvm::Div(momentumEqn_, rho_, U_, Uf_, options_.divScheme);                 // 对流项
                fvm::Laplacian(momentumEqn_, mu_, U_);                   // 扩散项
                fvm::Source(momentumEqn_, negGradPCell);           // 源项

                // 求解器设置
                Solver<Vector<Scalar>> momentumSolver(momentumEqn_,
                                                      options_.momentumMethod,
                                                      options_.momentumMaxIterations);

                momentumSolver.init(U_.getCellField().getData());
                if (options_.alphaU < 1.0 - 1e-6) {
                    momentumSolver.relax(options_.alphaU); // 松弛因子
                }
                momentumSolver.setTolerance(options_.momentumTolerance);
                momentumSolver.solve(); // 求解
                // 判断相对误差
                if (computeRelativeChange(
                    U_.getCellField().getData(),
                    U_.getCellField_0().getData()) <
                    options_.relativeTolerance) {
                    break;
                }
            }
        }
        inline void SIMPLE::solverPressureCorrEqn() {
            pCorr_.setValue(options_.pressureReferenceValue);
            // 计算负的速度散度
            CellField<Scalar> negDivUCell("-divU", mesh_, Scalar());
            std::vector<Scalar>& negDivUCellData = negDivUCell.getData();
            negDivUCellData = divergence(Uf_).getData();
            for (Scalar& divU : negDivUCellData) {
                divU = -divU;
            }
            for (int i = 0; i < 50; ++i) {
                pressureCorrEqn_.clear();
                // 扩散项, 源项
                fvm::Laplacian(pressureCorrEqn_, rAU_.getFaceField(), pCorr_);
                fvm::Source(pressureCorrEqn_, negDivUCell);
                // 若没有压力出口则使用参考点的压力，施加惩罚项
                if (options_.fixedPressurePatches.empty()) {
                    const ULL refCell = options_.pressureReferenceCell;
                    const Scalar penalty = options_.pressureReferencePenalty;
                    pressureCorrEqn_(refCell, refCell) += penalty;
                    pressureCorrEqn_.addB(refCell, penalty * options_.pressureReferenceValue);
                }
                Solver<Scalar> pressureCorrSolver(pressureCorrEqn_,
                                                  options_.pressureMethod,
                                                  options_.pressureMaxIterations);
                pressureCorrSolver.init(pCorr_.getCellField().getData());
                pressureCorrSolver.setTolerance(options_.pressureTolerance);
                pressureCorrSolver.solve();
                if (computeRelativeChange(
                    pCorr_.getCellField().getData(),
                    pCorr_.getCellField_0().getData()) <
                    options_.relativeTolerance ) {
                    break;
                }
            }
        }
        inline void SIMPLE::computeRAU() {
            std::vector<Scalar>& rAUCellFieldData = rAU_.getCellField().getData();
            const std::vector<Cell>& cells = mesh_->getCells();

            for (ULL cellId = 0; cellId < mesh_->getCellNumber(); ++cellId) {
                rAUCellFieldData[cellId] = cells[cellId].getVolume() / momentumEqn_(cellId, cellId);
            }
            rAU_.cellToFace();
        }
        inline void SIMPLE::initializepCorr() {
            pCorr_.initialize();
            initializeFixedPressureFaces();
            pCorr_.cellToFace();
        }
        inline void SIMPLE::initializeFixedPressureFaces() {
            std::unordered_set<std::string> patchNames(
                options_.fixedPressurePatches.begin(),
                options_.fixedPressurePatches.end()
            );

            const auto& patches = mesh_->getBoundaryPatches();

            // 确保所有给定的 patch 名称都存在于网格边界中
            if (!options_.fixedPressurePatches.empty()) {
                for (const std::string& name : patchNames) {
                    if (patches.find(name) == patches.end()) {
                        throw std::runtime_error(
                            "SIMPLE Error: fixed pressure patch \"" +
                            name +
                            "\" does not exist."
                        );
                    }
                }
            }

            for (const auto& [name, patch] : patches) {
                if (patch.getType() == BoundaryPatch::BoundaryType::EMPTY) {
                    continue;
                }
                if (patchNames.find(name) == patchNames.end()) {
                    pCorr_.setBoundaryCondition(name, 0, 1, 0); // 设置为0梯度边界条件
                    continue;
                }
                else {
                    pCorr_.setBoundaryCondition(name, 1, 0, 0); // 设置为固定压力边界条件;
                }
                const ULL startFace = patch.getStartFace();
                const ULL nFace = patch.getNFace();
                for (ULL faceId = startFace;
                     faceId < startFace + nFace;
                     ++faceId) {
                    fixedPressureFace_[faceId] = true;
                }
            }
        }
        inline void SIMPLE::correctFaceVelocity() {
            const std::vector<Face>& faces = mesh_->getFaces();
            const std::vector<Cell>& cells = mesh_->getCells();
            const CellField<Scalar>& pCorrCell = pCorr_.getCellField();
            const CellField<Vector<Scalar>>& gradPCorrCell = pCorr_.getCellGradientField();
            const CellField<Scalar>& rAUCell = rAU_.getCellField();

            for (ULL faceId : mesh_->getInternalFaceIndexes())
            {
                const Face& face = faces[faceId];
                const ULL owner =face.getOwnerIndex();
                const LL neighbor =face.getNeighborIndex();

                const Vector<Scalar> Sf =face.getArea() * face.getNormal();
                const Scalar magSf2 =Sf.magnitudeSquared();
                const Vector<Scalar> dPN =
                    cells[static_cast<ULL>(neighbor)].getCenter() -
                    cells[owner].getCenter();
                const Scalar dPNMag =dPN.magnitude();
                const Vector<Scalar> ePN =dPN / dPNMag;
                const Scalar denom = std::abs(face.getNormal() & ePN);
                const Scalar EfMag =face.getArea() / denom;
                const Scalar alpha =interpolationAlphaInternal(faceId);

                const Scalar Df =
                    util::interpolate(
                        rAUCell[owner],
                        rAUCell[static_cast<ULL>(neighbor)],
                        interpolation::Scheme::LINEAR,
                        alpha
                    );

                Scalar volumeFluxCorrection =
                    Df * EfMag / dPNMag *
                    (
                        pCorrCell[owner] -
                        pCorrCell[static_cast<ULL>(neighbor)]
                        );

                Uf_[faceId] +=
                    (volumeFluxCorrection / magSf2) * Sf;
            }
        }
        inline void SIMPLE::correctPressure() {
            std::vector<Scalar>& pData = p_.getCellField().getData();
            const std::vector<Scalar>& pCorrData = pCorr_.getCellField().getData();
            for (ULL cellId = 0; cellId < mesh_->getCellNumber(); ++cellId) {
                pData[cellId] += options_.alphaP * pCorrData[cellId];
            }
            p_.cellToFace();
        }
        inline void SIMPLE::correctVelocity() {
            std::vector<Vector<Scalar>>& UData = U_.getCellField().getData();
            const std::vector<Scalar>& rAUData = rAU_.getCellField().getData();
            const std::vector<Vector<Scalar>>& gardPCorrData = pCorr_.getCellGradientField().getData();
            for (ULL cellId = 0; cellId < mesh_->getCellNumber(); ++cellId) {
                UData[cellId] -= rAUData[cellId] * gardPCorrData[cellId];
            }
            U_.cellToFace();
        }
        inline Scalar SIMPLE::computeContinuityResidual() const {
            std::vector<Scalar> cellMassResidual(mesh_->getCellNumber(), Scalar());
            const std::vector<Face>& faces = mesh_->getFaces();

            Scalar totalAbsFlux = 0.0;
            // const std::vector<Vector<Scalar>>& UfData = U_.getFaceField().getData();
            const std::vector<Vector<Scalar>>& UfData = Uf_.getData();

            // 内部面遍历
            for (const ULL faceId : mesh_->getInternalFaceIndexes()) {
                const Face& face = faces[faceId];
                const ULL owner = face.getOwnerIndex();
                const LL neighbor = face.getNeighborIndex();
                const Scalar fluxValue = UfData[faceId] & face.getNormal() * face.getArea();

                cellMassResidual[owner] += fluxValue;
                cellMassResidual[neighbor] -= fluxValue;
                totalAbsFlux += std::abs(fluxValue);
            }
            // 边界面
            for (const ULL faceId : mesh_->getBoundaryFaceIndexes()) {
                const Face& face = faces[faceId];
                const ULL owner = face.getOwnerIndex();
                const Scalar fluxValue = UfData[faceId] & face.getNormal() * face.getArea();

                cellMassResidual[owner] += fluxValue;
                totalAbsFlux += std::abs(fluxValue);
            }
            // 找到cellMassResidual里最大值
            Scalar maxResidual = 0.0;
            for (const Scalar& residual : cellMassResidual) {
                maxResidual = std::max(maxResidual, std::abs(residual));
            }

            const Scalar avgAbsFlux = totalAbsFlux / mesh_->getCellNumber();

            return maxResidual / (avgAbsFlux + SMALL);
        }
        inline Scalar SIMPLE::computeRelativeChange(const std::vector<Scalar>& current, const std::vector<Scalar>& old) const {
            Scalar numerator = 0.0;
            Scalar denominator = 0.0;

            for (ULL i = 0; i < current.size(); ++i) {
                const Scalar diff = current[i] - old[i];

                numerator += diff * diff;
                denominator += current[i] * current[i];
            }

            return std::sqrt(numerator / (denominator + SMALL));
        }
        inline Scalar SIMPLE::computeRelativeChange(const std::vector<Vector<Scalar>>& current, const std::vector<Vector<Scalar>>& old) const
        {
            Scalar numerator = 0.0;
            Scalar denominator = 0.0;
            for (ULL i = 0; i < current.size(); ++i) {
                numerator += (current[i] - old[i]).magnitudeSquared();
                denominator += current[i].magnitudeSquared();
            }
            return std::sqrt(numerator / (denominator + SMALL));
        }
    }
}
#endif // SIMPLE1_H_