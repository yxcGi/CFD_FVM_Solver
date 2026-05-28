#ifndef SIMPLE_ALGORITHM_SIMPLE_HPP_
#define SIMPLE_ALGORITHM_SIMPLE_HPP_

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

namespace simple {

    class SIMPLE {
    public:
        using ULL = unsigned long long;
        using LL = long long;

        struct Options {
            // SIMPLE 外迭代控制
            int maxOuterIterations = 500;
            Scalar convergenceTolerance = 1e-8;
            int logInterval = 20;

            // 动量方程内部迭代
            int momentumMaxIterations = 5000;
            Scalar momentumTolerance = 1e-10;
            Solver<Vector<Scalar>>::Method momentumMethod =
                Solver<Vector<Scalar>>::Method::Jacobi;

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

            // 动量方程压力梯度源项系数。
            // 如果 p 是运动压力 p/rho，通常取 1。
            // 如果你传入的是物理压力，同时动量方程用的是 rho*U*U 和 mu，则也可取 1。
            Scalar pressureGradientCoefficient = 1.0;

            // 是否使用 Rhie-Chow 修正面通量
            bool useRhieChow = true;

            // 是否把修正后的质量通量投影回 U 的面场法向分量。
            // 这样下一次 fvm::Div() 用到的 FaceField<Vector> 与 SIMPLE 通量一致。
            bool projectFluxToFaceVelocity = true;

            // 压力修正边界 p'=0 的 patch 名称。
            // 比如压力出口可填 {"outlet"}。
            // 未列出的边界默认按 dp'/dn=0 处理。
            std::vector<std::string> fixedPressurePatches;
        };

        struct Residual {
            int iteration = 0;
            Scalar continuity = std::numeric_limits<Scalar>::max();
            Scalar velocity = std::numeric_limits<Scalar>::max();
            Scalar pressure = std::numeric_limits<Scalar>::max();
            Scalar pressureEquation = std::numeric_limits<Scalar>::max();

            bool converged(Scalar tolerance) const {
                return continuity < tolerance &&
                    velocity < tolerance &&
                    pressure < tolerance;
            }
        };

    public:
        SIMPLE(
            Field<Vector<Scalar>> &U,
            Field<Scalar> &p,
            FaceField<Scalar> &rho,
            FaceField<Scalar> &gamma
        );

        SIMPLE(
            Field<Vector<Scalar>> &U,
            Field<Scalar> &p,
            FaceField<Scalar> &rho,
            FaceField<Scalar> &gamma,
            const Options &options
        );

        Residual solve() {
            // 初始插值，保证 U、p 的面场和梯度场有效。
            U_.cellToFace();
            p_.cellToFace();

            computeRhieChowFlux();

            if (options_.projectFluxToFaceVelocity) {
                projectFluxToUFaceField();
            }

            Residual residual;

            for (int outerIter = 0;
                 outerIter < options_.maxOuterIterations;
                 ++outerIter) {
                residual.iteration = outerIter + 1;

                const std::vector<Vector<Scalar>> oldU =
                    U_.getCellField().getData();
                const std::vector<Scalar> oldP =
                    p_.getCellField().getData();

                solveMomentumEquation();

                // 动量方程求出 U* 后，重新计算 U* 面场和 p 梯度。
                U_.cellToFace();
                p_.cellToFace();

                computeRhieChowFlux();

                residual.pressureEquation = solvePressureCorrectionEquation();

                const std::vector<Vector<Scalar>> gradPCorr =
                    computeScalarCellGradient(pCorr_.getCellField().getData());

                correctPressure();
                correctVelocity(gradPCorr);

                // pressure correction equation 已经修正了 phi_。
                // 这里更新 U、p 的面场和梯度，再把 phi_ 投影回 U 面法向分量。
                U_.cellToFace();
                p_.cellToFace();

                if (options_.projectFluxToFaceVelocity) {
                    projectFluxToUFaceField();
                }

                residual.continuity = computeContinuityResidual();
                residual.velocity = computeRelativeChange(
                    U_.getCellField().getData(),
                    oldU
                );
                residual.pressure = computeRelativeChange(
                    p_.getCellField().getData(),
                    oldP
                );

                if (options_.logInterval > 0 &&
                    ((outerIter + 1) % options_.logInterval == 0 ||
                     outerIter == 0 ||
                     residual.converged(options_.convergenceTolerance))) {
                    std::cout
                        << "SIMPLE outer iteration " << outerIter + 1
                        << ", continuity = " << residual.continuity
                        << ", U = " << residual.velocity
                        << ", p = " << residual.pressure
                        << ", pEqn = " << residual.pressureEquation
                        << std::endl;
                }

                if (residual.converged(options_.convergenceTolerance)) {
                    break;
                }
            }

            return residual;
        }

        const FaceField<Scalar> &flux() const {
            return phi_;
        }

        FaceField<Scalar> &flux() {
            return phi_;
        }

        const Field<Scalar> &pressureCorrection() const {
            return pCorr_;
        }

        const std::vector<Scalar> &rAU() const {
            return rAU_;
        }

    private:
        void validateInput() const {
            if (mesh_ == nullptr) {
                throw std::runtime_error("SIMPLE Error: mesh pointer is null.");
            }

            if (!mesh_->isValid()) {
                throw std::runtime_error("SIMPLE Error: mesh is not valid.");
            }

            if (p_.getMesh() != mesh_ ||
                rho_.getMesh() != mesh_ ||
                gamma_.getMesh() != mesh_) {
                throw std::runtime_error(
                    "SIMPLE Error: U, p, rho, gamma must use the same mesh."
                );
            }

            if (!U_.isValid()) {
                throw std::runtime_error("SIMPLE Error: U field is not valid.");
            }

            if (!p_.isValid()) {
                throw std::runtime_error("SIMPLE Error: p field is not valid.");
            }

            if (!rho_.isValid()) {
                throw std::runtime_error("SIMPLE Error: rho face field is not valid.");
            }

            if (!gamma_.isValid()) {
                throw std::runtime_error("SIMPLE Error: gamma face field is not valid.");
            }
        }

        void validateOptions() const {
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

        void initializeFixedPressureFaces() {
            if (options_.fixedPressurePatches.empty()) {
                return;
            }

            std::unordered_set<std::string> patchNames(
                options_.fixedPressurePatches.begin(),
                options_.fixedPressurePatches.end()
            );

            const auto &patches = mesh_->getBoundaryPatches();

            for (const std::string &name : patchNames) {
                if (patches.find(name) == patches.end()) {
                    throw std::runtime_error(
                        "SIMPLE Error: fixed pressure patch \"" +
                        name +
                        "\" does not exist."
                    );
                }
            }

            for (const auto &[name, patch] : patches) {
                if (patchNames.find(name) == patchNames.end()) {
                    continue;
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

        void solveMomentumEquation() {
            momentumEqn_.clear();

            // 对流项：进入矩阵左端。
            fvm::Div(
                momentumEqn_,
                rho_,
                U_,
                U_.getFaceField(),
                options_.divScheme
            );

            // 扩散项：按照你现有实现，进入矩阵左端。
            // 这里 gamma 可以理解为 mu 或 nu，取决于你整体方程的量纲约定。
            fvm::Laplacian(momentumEqn_, gamma_, U_);

            addPressureGradientSource(momentumEqn_);

            const std::vector<Vector<Scalar>> oldU =
                U_.getCellField().getData();

            if (options_.alphaU < 1.0) {
                relaxMomentumEquation(momentumEqn_, oldU, options_.alphaU);
            }

            computeRAU(momentumEqn_);

            Solver<Vector<Scalar>> solver(
                momentumEqn_,
                options_.momentumMethod,
                options_.momentumMaxIterations
            );

            solver.init(U_.getCellField().getData());
            solver.setTolerance(options_.momentumTolerance);
            solver.solve();
        }

        void addPressureGradientSource(
            SparseMatrix<Vector<Scalar>> &equation
        ) {
            const CellField<Vector<Scalar>> &gradP =
                p_.getCellGradientField();

            const std::vector<Cell> &cells = mesh_->getCells();

            for (ULL cellId = 0;
                 cellId < mesh_->getCellNumber();
                 ++cellId) {
                const Vector<Scalar> pressureSource =
                    -options_.pressureGradientCoefficient *
                    cells[cellId].getVolume() *
                    gradP[cellId];

                equation.addB(cellId, pressureSource);
            }
        }

        void relaxMomentumEquation(
            SparseMatrix<Vector<Scalar>> &equation,
            const std::vector<Vector<Scalar>> &oldU,
            Scalar alpha
        ) {
            for (ULL cellId = 0;
                 cellId < mesh_->getCellNumber();
                 ++cellId) {
                const Scalar oldDiag = equation(cellId, cellId);

                if (std::abs(oldDiag) < SMALL) {
                    throw std::runtime_error(
                        "SIMPLE Error: zero diagonal in momentum equation."
                    );
                }

                const Scalar relaxedDiag = oldDiag / alpha;

                equation.setValue(cellId, cellId, relaxedDiag);

                equation.addB(
                    cellId,
                    (1.0 - alpha) * relaxedDiag * oldU[cellId]
                );
            }
        }

        void computeRAU(const SparseMatrix<Vector<Scalar>> &equation) {
            const std::vector<Cell> &cells = mesh_->getCells();

            for (ULL cellId = 0;
                 cellId < mesh_->getCellNumber();
                 ++cellId) {
                const Scalar diag = equation.at(cellId, cellId);

                if (std::abs(diag) < SMALL) {
                    throw std::runtime_error(
                        "SIMPLE Error: zero diagonal when computing rAU."
                    );
                }

                // 这里 rAU_ 实际保存 D_P = V_P / a_P。
                rAU_[cellId] = cells[cellId].getVolume() / diag;
            }
        }

        Scalar solvePressureCorrectionEquation() {
            std::vector<Scalar> &pCorr =
                pCorr_.getCellField().getData();

            std::fill(pCorr.begin(), pCorr.end(), Scalar{});
            pCorr_.getCellField_0().getData() = pCorr;

            Scalar pressureEquationResidual =
                std::numeric_limits<Scalar>::max();

            std::vector<Vector<Scalar>> gradPCorr(
                mesh_->getCellNumber(),
                Vector<Scalar>{}
            );

            for (int nonOrth = 0;
                 nonOrth < options_.nNonOrthogonalCorrectors;
                 ++nonOrth) {
                const bool useExplicitNonOrth =
                    options_.nNonOrthogonalCorrectors > 1;

                assemblePressureCorrectionEquation(
                    gradPCorr,
                    useExplicitNonOrth
                );

                pressureEquationResidual = solveScalarEquationJacobi(
                    pressureCorrEqn_,
                    pCorr,
                    options_.pressureMaxIterations,
                    options_.pressureTolerance
                );

                pCorr_.getCellField().getData() = pCorr;
                pCorr_.getCellField_0().getData() = pCorr;

                gradPCorr = computeScalarCellGradient(pCorr);

                if (nonOrth == options_.nNonOrthogonalCorrectors - 1) {
                    correctFlux(gradPCorr, useExplicitNonOrth);
                }
            }

            return pressureEquationResidual;
        }

        void assemblePressureCorrectionEquation(
            const std::vector<Vector<Scalar>> &gradPCorr,
            bool useExplicitNonOrth
        ) {
            pressureCorrEqn_.clear();

            const std::vector<Face> &faces = mesh_->getFaces();
            const std::vector<ULL> &internalFaceIndexes =
                mesh_->getInternalFaceIndexes();
            const std::vector<ULL> &boundaryFaceIndexes =
                mesh_->getBoundaryFaceIndexes();

            for (ULL faceId : internalFaceIndexes) {
                if (isEmptyFace(faceId)) {
                    continue;
                }

                const Face &face = faces[faceId];

                const ULL owner = face.getOwnerIndex();
                const LL neighbor = face.getNeighborIndex();

                if (neighbor < 0) {
                    continue;
                }

                const Scalar coeff = pressureCorrectionCoeffInternal(faceId);

                pressureCorrEqn_.addValue(owner, owner, coeff);
                pressureCorrEqn_.addValue(owner, static_cast<ULL>(neighbor), -coeff);

                pressureCorrEqn_.addValue(static_cast<ULL>(neighbor),
                                          static_cast<ULL>(neighbor),
                                          coeff);
                pressureCorrEqn_.addValue(static_cast<ULL>(neighbor), owner, -coeff);

                Scalar explicitNonOrthFlux = 0.0;

                if (useExplicitNonOrth) {
                    explicitNonOrthFlux =
                        pressureCorrectionNonOrthogonalFluxInternal(
                            faceId,
                            gradPCorr
                        );
                }

                // 连续性：sum(phi + phi') = 0。
                // 对 owner，phi 正方向为 owner -> neighbor。
                const Scalar rhsFlux = phi_[faceId] + explicitNonOrthFlux;

                pressureCorrEqn_.addB(owner, -rhsFlux);
                pressureCorrEqn_.addB(static_cast<ULL>(neighbor), rhsFlux);
            }

            for (ULL faceId : boundaryFaceIndexes) {
                if (isEmptyFace(faceId)) {
                    continue;
                }

                const Face &face = faces[faceId];
                const ULL owner = face.getOwnerIndex();

                if (fixedPressureFace_[faceId]) {
                    const Scalar coeff =
                        pressureCorrectionCoeffBoundary(faceId);

                    pressureCorrEqn_.addValue(owner, owner, coeff);
                }

                // 对 dp'/dn = 0 的边界，只修正内部通量，不修正边界通量。
                // 对 p'=0 的压力边界，上面已加入边界压力修正系数。
                pressureCorrEqn_.addB(owner, -phi_[faceId]);
            }

            applyPressureReference(pressureCorrEqn_);
        }

        void applyPressureReference(SparseMatrix<Scalar> &equation) {
            const ULL refCell = options_.pressureReferenceCell;
            const Scalar penalty = options_.pressureReferencePenalty;

            equation.addValue(refCell, refCell, penalty);
            equation.addB(
                refCell,
                penalty * options_.pressureReferenceValue
            );
        }

        Scalar solveScalarEquationJacobi(
            const SparseMatrix<Scalar> &equation,
            std::vector<Scalar> &x,
            int maxIterations,
            Scalar tolerance
        ) const {
            const ULL size = equation.size();

            if (x.size() != size) {
                x.assign(size, Scalar{});
            }

            const std::vector<Scalar> &values = equation.getValues();
            const std::vector<ULL> &colIndex = equation.getColIndexs();
            const std::vector<ULL> &rowPtr = equation.getRowPointer();
            const std::vector<Scalar> &b = equation.getB();

            std::vector<Scalar> xOld = x;
            std::vector<Scalar> xNew = x;

            std::vector<Scalar> diag(size, Scalar{});

            for (ULL row = 0; row < size; ++row) {
                diag[row] = equation.at(row, row);

                if (std::abs(diag[row]) < SMALL) {
                    throw std::runtime_error(
                        "SIMPLE Error: zero diagonal in scalar Jacobi solver."
                    );
                }
            }

            Scalar maxResidual = std::numeric_limits<Scalar>::max();

            for (int iter = 0; iter < maxIterations; ++iter) {
                for (ULL row = 0; row < size; ++row) {
                    Scalar sum = 0.0;

                    for (ULL id = rowPtr[row];
                         id < rowPtr[row + 1];
                         ++id) {
                        const ULL col = colIndex[id];

                        if (col == row) {
                            continue;
                        }

                        sum += values[id] * xOld[col];
                    }

                    xNew[row] = (b[row] - sum) / diag[row];
                }

                maxResidual = 0.0;

                for (ULL row = 0; row < size; ++row) {
                    Scalar Ax = 0.0;

                    for (ULL id = rowPtr[row];
                         id < rowPtr[row + 1];
                         ++id) {
                        Ax += values[id] * xNew[colIndex[id]];
                    }

                    maxResidual =
                        std::max(maxResidual, std::abs(b[row] - Ax));
                }

                xOld.swap(xNew);

                if (maxResidual < tolerance) {
                    break;
                }
            }

            x = xOld;

            return maxResidual;
        }

        void correctPressure() {
            CellField<Scalar> &pCell = p_.getCellField();
            CellField<Scalar> &pOld = p_.getCellField_0();
            const CellField<Scalar> &pCorrCell = pCorr_.getCellField();

            for (ULL cellId = 0;
                 cellId < mesh_->getCellNumber();
                 ++cellId) {
                pCell[cellId] += options_.alphaP * pCorrCell[cellId];
            }

            pOld.getData() = pCell.getData();
        }

        void correctVelocity(
            const std::vector<Vector<Scalar>> &gradPCorr
        ) {
            CellField<Vector<Scalar>> &UCell = U_.getCellField();
            CellField<Vector<Scalar>> &UOld = U_.getCellField_0();

            for (ULL cellId = 0;
                 cellId < mesh_->getCellNumber();
                 ++cellId) {
                UCell[cellId] -= rAU_[cellId] * gradPCorr[cellId];
            }

            UOld.getData() = UCell.getData();
        }

        void correctFlux(
            const std::vector<Vector<Scalar>> &gradPCorr,
            bool useExplicitNonOrth
        ) {
            const std::vector<Face> &faces = mesh_->getFaces();
            const std::vector<ULL> &internalFaceIndexes =
                mesh_->getInternalFaceIndexes();
            const std::vector<ULL> &boundaryFaceIndexes =
                mesh_->getBoundaryFaceIndexes();

            const CellField<Scalar> &pCorrCell =
                pCorr_.getCellField();

            for (ULL faceId : internalFaceIndexes) {
                if (isEmptyFace(faceId)) {
                    phi_.setValue(faceId, Scalar{});
                    continue;
                }

                const Face &face = faces[faceId];

                const ULL owner = face.getOwnerIndex();
                const LL neighbor = face.getNeighborIndex();

                if (neighbor < 0) {
                    continue;
                }

                const Scalar coeff =
                    pressureCorrectionCoeffInternal(faceId);

                Scalar fluxCorrection =
                    coeff *
                    (pCorrCell[owner] -
                     pCorrCell[static_cast<ULL>(neighbor)]);

                if (useExplicitNonOrth) {
                    fluxCorrection +=
                        pressureCorrectionNonOrthogonalFluxInternal(
                            faceId,
                            gradPCorr
                        );
                }

                phi_.setValue(faceId, phi_[faceId] + fluxCorrection);
            }

            for (ULL faceId : boundaryFaceIndexes) {
                if (isEmptyFace(faceId)) {
                    phi_.setValue(faceId, Scalar{});
                    continue;
                }

                if (!fixedPressureFace_[faceId]) {
                    continue;
                }

                const Face &face = faces[faceId];
                const ULL owner = face.getOwnerIndex();

                const Scalar coeff =
                    pressureCorrectionCoeffBoundary(faceId);

                phi_.setValue(
                    faceId,
                    phi_[faceId] + coeff * pCorrCell[owner]
                );
            }
        }

        void computeRhieChowFlux() {
            const std::vector<Face> &faces = mesh_->getFaces();
            const std::vector<ULL> &internalFaceIndexes =
                mesh_->getInternalFaceIndexes();
            const std::vector<ULL> &boundaryFaceIndexes =
                mesh_->getBoundaryFaceIndexes();

            const CellField<Vector<Scalar>> &UCell = U_.getCellField();
            const FaceField<Vector<Scalar>> &UFace = U_.getFaceField();

            const CellField<Scalar> &pCell = p_.getCellField();
            const CellField<Vector<Scalar>> &gradP =
                p_.getCellGradientField();

            for (ULL faceId : internalFaceIndexes) {
                if (isEmptyFace(faceId)) {
                    phi_.setValue(faceId, Scalar{});
                    continue;
                }

                const Face &face = faces[faceId];

                const ULL owner = face.getOwnerIndex();
                const LL neighbor = face.getNeighborIndex();

                if (neighbor < 0) {
                    continue;
                }

                const Vector<Scalar> Sf =
                    face.getArea() * face.getNormal();

                if (!options_.useRhieChow) {
                    phi_.setValue(
                        faceId,
                        rho_[faceId] * (UFace[faceId] & Sf)
                    );
                    continue;
                }

                const Scalar alpha = interpolationAlphaInternal(faceId);

                const Vector<Scalar> Uf =
                    (1.0 - alpha) * UCell[owner] +
                    alpha * UCell[static_cast<ULL>(neighbor)];

                const Vector<Scalar> gradPF =
                    (1.0 - alpha) * gradP[owner] +
                    alpha * gradP[static_cast<ULL>(neighbor)];

                const Scalar Df =
                    (1.0 - alpha) * rAU_[owner] +
                    alpha * rAU_[static_cast<ULL>(neighbor)];

                const Scalar gradPFaceDotSf =
                    faceNormalPressureGradientDotSf(
                        faceId,
                        pCell[owner],
                        pCell[static_cast<ULL>(neighbor)],
                        gradPF
                    );

                const Scalar massFlux =
                    rho_[faceId] *
                    (
                        (Uf & Sf) +
                        Df * ((gradPF & Sf) - gradPFaceDotSf)
                        );

                phi_.setValue(faceId, massFlux);
            }

            for (ULL faceId : boundaryFaceIndexes) {
                if (isEmptyFace(faceId)) {
                    phi_.setValue(faceId, Scalar{});
                    continue;
                }

                const Face &face = faces[faceId];
                const Vector<Scalar> Sf =
                    face.getArea() * face.getNormal();

                phi_.setValue(
                    faceId,
                    rho_[faceId] * (UFace[faceId] & Sf)
                );
            }
        }

        void projectFluxToUFaceField() {
            FaceField<Vector<Scalar>> &UFace = U_.getFaceField();
            const std::vector<Face> &faces = mesh_->getFaces();

            for (ULL faceId = 0;
                 faceId < mesh_->getFaceNumber();
                 ++faceId) {
                if (isEmptyFace(faceId)) {
                    UFace.setValue(faceId, Vector<Scalar>{});
                    continue;
                }

                const Scalar rhoF = rho_[faceId];

                if (std::abs(rhoF) < SMALL) {
                    continue;
                }

                const Face &face = faces[faceId];
                const Vector<Scalar> normal = face.getNormal();
                const Vector<Scalar> Sf = face.getArea() * normal;

                const Scalar targetVolumeFlux = phi_[faceId] / rhoF;
                const Scalar currentVolumeFlux = UFace[faceId] & Sf;

                const Scalar normalVelocityCorrection =
                    (targetVolumeFlux - currentVolumeFlux) /
                    std::max(face.getArea(), SMALL);

                UFace.setValue(
                    faceId,
                    UFace[faceId] + normalVelocityCorrection * normal
                );
            }
        }

        std::vector<Vector<Scalar>> computeScalarCellGradient(
            const std::vector<Scalar> &scalarCellValue
        ) const {
            if (scalarCellValue.size() != mesh_->getCellNumber()) {
                throw std::runtime_error(
                    "SIMPLE Error: scalarCellValue size is inconsistent with mesh."
                );
            }

            const std::vector<Cell> &cells = mesh_->getCells();
            const std::vector<Face> &faces = mesh_->getFaces();

            std::vector<Vector<Scalar>> grad(
                mesh_->getCellNumber(),
                Vector<Scalar>{}
            );

            for (ULL cellId = 0;
                 cellId < mesh_->getCellNumber();
                 ++cellId) {
                const Cell &cell = cells[cellId];

                Vector<Scalar> total{};

                for (ULL faceId : cell.getFaceIndexes()) {
                    if (isEmptyFace(faceId)) {
                        continue;
                    }

                    const Face &face = faces[faceId];

                    Scalar faceValue = 0.0;

                    const ULL owner = face.getOwnerIndex();
                    const LL neighbor = face.getNeighborIndex();

                    if (neighbor >= 0) {
                        const Scalar alpha =
                            interpolationAlphaInternal(faceId);

                        faceValue =
                            (1.0 - alpha) * scalarCellValue[owner] +
                            alpha * scalarCellValue[static_cast<ULL>(neighbor)];
                    }
                    else {
                        if (fixedPressureFace_[faceId]) {
                            faceValue = 0.0;
                        }
                        else {
                            faceValue = scalarCellValue[owner];
                        }
                    }

                    const Vector<Scalar> Sf =
                        face.getArea() * face.getNormal();

                    if (owner == cellId) {
                        total += faceValue * Sf;
                    }
                    else {
                        total -= faceValue * Sf;
                    }
                }

                grad[cellId] = total / cell.getVolume();
            }

            return grad;
        }

        Scalar pressureCorrectionCoeffInternal(ULL faceId) const {
            const std::vector<Face> &faces = mesh_->getFaces();
            const std::vector<Cell> &cells = mesh_->getCells();

            const Face &face = faces[faceId];

            const ULL owner = face.getOwnerIndex();
            const LL neighbor = face.getNeighborIndex();

            if (neighbor < 0) {
                throw std::runtime_error(
                    "SIMPLE Error: internal pressure coefficient requested on boundary face."
                );
            }

            const Vector<Scalar> d =
                cells[static_cast<ULL>(neighbor)].getCenter() -
                cells[owner].getCenter();

            const Scalar distance = d.magnitude();

            if (distance < SMALL) {
                throw std::runtime_error(
                    "SIMPLE Error: zero owner-neighbor distance."
                );
            }

            const Vector<Scalar> e = d / distance;

            const Scalar denom = std::abs(face.getNormal() & e);

            if (denom < SMALL) {
                throw std::runtime_error(
                    "SIMPLE Error: bad non-orthogonal face geometry."
                );
            }

            const Scalar EfMag = face.getArea() / denom;

            const Scalar alpha = interpolationAlphaInternal(faceId);

            const Scalar Df =
                (1.0 - alpha) * rAU_[owner] +
                alpha * rAU_[static_cast<ULL>(neighbor)];

            return rho_[faceId] * Df * EfMag / distance;
        }

        Scalar pressureCorrectionCoeffBoundary(ULL faceId) const {
            const std::vector<Face> &faces = mesh_->getFaces();
            const std::vector<Cell> &cells = mesh_->getCells();

            const Face &face = faces[faceId];

            const ULL owner = face.getOwnerIndex();

            const Vector<Scalar> d =
                face.getCenter() -
                cells[owner].getCenter();

            const Scalar distance = d.magnitude();

            if (distance < SMALL) {
                throw std::runtime_error(
                    "SIMPLE Error: zero cell-boundary distance."
                );
            }

            const Vector<Scalar> e = d / distance;

            const Scalar denom = std::abs(face.getNormal() & e);

            if (denom < SMALL) {
                throw std::runtime_error(
                    "SIMPLE Error: bad boundary face geometry."
                );
            }

            const Scalar EfMag = face.getArea() / denom;

            return rho_[faceId] * rAU_[owner] * EfMag / distance;
        }

        Scalar pressureCorrectionNonOrthogonalFluxInternal(
            ULL faceId,
            const std::vector<Vector<Scalar>> &gradPCorr
        ) const {
            const std::vector<Face> &faces = mesh_->getFaces();
            const std::vector<Cell> &cells = mesh_->getCells();

            const Face &face = faces[faceId];

            const ULL owner = face.getOwnerIndex();
            const LL neighbor = face.getNeighborIndex();

            if (neighbor < 0) {
                return 0.0;
            }

            const Scalar alpha = interpolationAlphaInternal(faceId);

            const Vector<Scalar> gradPF =
                (1.0 - alpha) * gradPCorr[owner] +
                alpha * gradPCorr[static_cast<ULL>(neighbor)];

            const Scalar Df =
                (1.0 - alpha) * rAU_[owner] +
                alpha * rAU_[static_cast<ULL>(neighbor)];

            const Vector<Scalar> d =
                cells[static_cast<ULL>(neighbor)].getCenter() -
                cells[owner].getCenter();

            const Vector<Scalar> e = d.unitVector();

            const Scalar denom = std::abs(face.getNormal() & e);

            if (denom < SMALL) {
                return 0.0;
            }

            const Vector<Scalar> Sf =
                face.getArea() * face.getNormal();

            const Scalar EfMag = face.getArea() / denom;
            const Vector<Scalar> Ef = EfMag * e;
            const Vector<Scalar> Tf = Sf - Ef;

            return -rho_[faceId] * Df * (gradPF & Tf);
        }

        Scalar faceNormalPressureGradientDotSf(
            ULL faceId,
            Scalar pOwner,
            Scalar pNeighbor,
            const Vector<Scalar> &gradPF
        ) const {
            const std::vector<Face> &faces = mesh_->getFaces();
            const std::vector<Cell> &cells = mesh_->getCells();

            const Face &face = faces[faceId];

            const ULL owner = face.getOwnerIndex();
            const LL neighbor = face.getNeighborIndex();

            if (neighbor < 0) {
                return gradPF & (face.getArea() * face.getNormal());
            }

            const Vector<Scalar> d =
                cells[static_cast<ULL>(neighbor)].getCenter() -
                cells[owner].getCenter();

            const Scalar distance = d.magnitude();

            if (distance < SMALL) {
                throw std::runtime_error(
                    "SIMPLE Error: zero owner-neighbor distance."
                );
            }

            const Vector<Scalar> e = d / distance;

            const Scalar denom = std::abs(face.getNormal() & e);

            if (denom < SMALL) {
                throw std::runtime_error(
                    "SIMPLE Error: bad face geometry in Rhie-Chow interpolation."
                );
            }

            const Vector<Scalar> Sf =
                face.getArea() * face.getNormal();

            const Scalar EfMag = face.getArea() / denom;
            const Vector<Scalar> Ef = EfMag * e;
            const Vector<Scalar> Tf = Sf - Ef;

            return ((pNeighbor - pOwner) / distance) * EfMag +
                (gradPF & Tf);
        }

        Scalar interpolationAlphaInternal(ULL faceId) const {
            const std::vector<Face> &faces = mesh_->getFaces();
            const std::vector<Cell> &cells = mesh_->getCells();

            const Face &face = faces[faceId];

            const ULL owner = face.getOwnerIndex();
            const LL neighbor = face.getNeighborIndex();

            if (neighbor < 0) {
                return 0.0;
            }

            const Vector<Scalar> &faceCenter = face.getCenter();
            const Vector<Scalar> &ownerCenter = cells[owner].getCenter();
            const Vector<Scalar> &neighborCenter =
                cells[static_cast<ULL>(neighbor)].getCenter();

            const Vector<Scalar> &normal = face.getNormal();

            const Scalar ownerDistance =
                (faceCenter - ownerCenter) & normal;
            const Scalar neighborDistance =
                (neighborCenter - faceCenter) & normal;

            const Scalar denominator =
                ownerDistance + neighborDistance;

            if (std::abs(denominator) < SMALL) {
                return 0.5;
            }

            Scalar alpha = ownerDistance / denominator;

            if (alpha < 0.0 || alpha > 1.0) {
                alpha = 0.5;
            }

            return alpha;
        }

        Scalar computeContinuityResidual() const {
            std::vector<Scalar> cellMassResidual(
                mesh_->getCellNumber(),
                Scalar{}
            );

            const std::vector<Face> &faces = mesh_->getFaces();

            Scalar totalAbsFlux = 0.0;

            for (ULL faceId = 0;
                 faceId < mesh_->getFaceNumber();
                 ++faceId) {
                if (isEmptyFace(faceId)) {
                    continue;
                }

                const Face &face = faces[faceId];

                const ULL owner = face.getOwnerIndex();
                const LL neighbor = face.getNeighborIndex();

                const Scalar fluxValue = phi_[faceId];

                cellMassResidual[owner] += fluxValue;

                if (neighbor >= 0) {
                    cellMassResidual[static_cast<ULL>(neighbor)] -= fluxValue;
                }

                totalAbsFlux += std::abs(fluxValue);
            }

            Scalar maxResidual = 0.0;

            for (Scalar r : cellMassResidual) {
                maxResidual = std::max(maxResidual, std::abs(r));
            }

            const Scalar fluxScale =
                totalAbsFlux /
                std::max<Scalar>(static_cast<Scalar>(mesh_->getCellNumber()), 1.0);

            return maxResidual / (fluxScale + SMALL);
        }

        Scalar computeRelativeChange(
            const std::vector<Scalar> &current,
            const std::vector<Scalar> &old
        ) const {
            if (current.size() != old.size()) {
                throw std::runtime_error(
                    "SIMPLE Error: scalar field size mismatch."
                );
            }

            Scalar numerator = 0.0;
            Scalar denominator = 0.0;

            for (ULL i = 0; i < current.size(); ++i) {
                const Scalar diff = current[i] - old[i];

                numerator += diff * diff;
                denominator += current[i] * current[i];
            }

            return std::sqrt(numerator / (denominator + SMALL));
        }

        Scalar computeRelativeChange(
            const std::vector<Vector<Scalar>> &current,
            const std::vector<Vector<Scalar>> &old
        ) const {
            if (current.size() != old.size()) {
                throw std::runtime_error(
                    "SIMPLE Error: vector field size mismatch."
                );
            }

            Scalar numerator = 0.0;
            Scalar denominator = 0.0;

            for (ULL i = 0; i < current.size(); ++i) {
                numerator += (current[i] - old[i]).magnitudeSquared();
                denominator += current[i].magnitudeSquared();
            }

            return std::sqrt(numerator / (denominator + SMALL));
        }

        bool isEmptyFace(ULL faceId) const {
            if (mesh_->getDimension() != Mesh::Dimension::TWO_D) {
                return false;
            }

            const auto &emptyPair = mesh_->getEmptyFaceIndexesPair();

            return faceId >= emptyPair.first &&
                faceId < emptyPair.second;
        }

    private:
        static constexpr Scalar SMALL = 1e-30;

    private:
        Mesh *mesh_;

        Field<Vector<Scalar>> &U_;
        Field<Scalar> &p_;

        FaceField<Scalar> &rho_;
        FaceField<Scalar> &gamma_;

        Options options_;

        FaceField<Scalar> phi_;
        Field<Scalar> pCorr_;

        SparseMatrix<Vector<Scalar>> momentumEqn_;
        SparseMatrix<Scalar> pressureCorrEqn_;

        // 实际保存 D_P = V_P / a_P。
        std::vector<Scalar> rAU_;

        std::vector<bool> fixedPressureFace_;
    };

    inline SIMPLE::SIMPLE(
        Field<Vector<Scalar>> &U,
        Field<Scalar> &p,
        FaceField<Scalar> &rho,
        FaceField<Scalar> &gamma
    )
        : SIMPLE(U, p, rho, gamma, Options{}) {}

    inline SIMPLE::SIMPLE(
        Field<Vector<Scalar>> &U,
        Field<Scalar> &p,
        FaceField<Scalar> &rho,
        FaceField<Scalar> &gamma,
        const Options &options
    )
        : mesh_(U.getMesh())
        , U_(U)
        , p_(p)
        , rho_(rho)
        , gamma_(gamma)
        , options_(options)
        , phi_("phi", mesh_, Scalar{})
        , pCorr_("pCorr", mesh_)
        , momentumEqn_(mesh_)
        , pressureCorrEqn_(mesh_) {
        validateInput();
        validateOptions();

        const ULL cellNumber = mesh_->getCellNumber();
        const ULL faceNumber = mesh_->getFaceNumber();

        rAU_.assign(cellNumber, Scalar{});
        fixedPressureFace_.assign(faceNumber, false);

        pCorr_.initialize();

        initializeFixedPressureFaces();
    }
    
} // namespace simple

#endif // SIMPLE_ALGORITHM_SIMPLE_HPP_