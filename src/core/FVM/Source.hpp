#ifndef SOURCE_H_
#define SOURCE_H_

#include "SparseMatrix.hpp"
#include "Field.hpp"
#include <cassert>
#include "interpolation.hpp"

namespace fvm {

    /// 源项：S_phi = Sp * phi + S_c
    
    /**
     * @brief 只含S_c的源项
     * @tparam Tp 场类型
     * @param matrix 待赋值矩阵
     * @param S_c 源项场
     */
    template <typename Tp>
    void Source(
        SparseMatrix<Tp> &matrix,
        const CellField<Tp> &S_c
    );

    template<typename Tp>
    void Source(SparseMatrix<Tp> &matrix,
                const CellField<Tp> &S_c)
    {
        Mesh *mesh = matrix.getMesh();

        // 判断参数是否有效
        if (!matrix.isValid()) {
            std::cerr << "Source() Error: matrix is not valid." << std::endl;
            throw std::invalid_argument("matrix is not valid.");
        }
        if (!S_c.isValid()) {
            std::cerr << "Source() Error: S_c is not valid." << std::endl;
            throw std::invalid_argument("Source is not valid.");
        }
        // 判断是否矩阵是否为空指针网格（稀疏矩阵只能由mesh构造）
        if (mesh == nullptr) {
            std::cerr << "Source() Error: matrix is not constructed by mesh." << std::endl;
            throw std::invalid_argument("matrix is not constructed by mesh.");
        }
        // 判断输入参数是否为同一网格
        if (mesh != S_c.getMesh()) {
            std::cerr << "Source() Error: matrix and S_c are not constructed by the same mesh." << std::endl;
            throw std::invalid_argument("matrix and S_c are not constructed by the same mesh.");
        }

        using ULL = unsigned long long;
        const std::vector<Cell> &cells = mesh->getCells();

        for (ULL cellId = 0; cellId < mesh->getCellNumber(); ++cellId) {
            matrix.addB(cellId, cells[cellId].getVolume() * S_c[cellId]);
        }
    }

}

#endif // SOURCE_H_