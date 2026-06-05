#ifndef SOLVER_H_
#define SOLVER_H_

#include "SparseMatrix.hpp"
#include "Field.hpp"
#include "threadpool.hpp"


template <typename Tp>
class Solver
{
    using ULL = unsigned long long;
    static constexpr double DEFAULT_EPSILON = 1e-6; // 默认的精度
    static constexpr int DEFAULT_OUTPUT_INTERVAL = 20;
public:
    // 求解方法
    enum class Method
    {
        Jacobi,             // 雅可比迭代
        GaussSeidel,        // 高斯赛德尔迭代
        AMG,                // 代数多重网格法
    };

public:
    Solver() = delete;
    explicit Solver(SparseMatrix<Tp>& matrix, Method method, int maxmaxIterationNum);
    Solver(const Solver&) = delete;
    Solver(Solver&&) = delete;
    ~Solver() = default;

public:
    // 初始化
    void init(const std::vector<Tp>& x0);
    void init(Tp value0);
    void init(const Field<Tp>& initField);
    // 获取初始化后的残差，判据：两次外迭代解的相对变化量
    Scalar Error();

    // 并行
    void setParallel();

    // 求解方程
    void solve();

    // 设置容差
    void setTolerance(Scalar tolerance);

    // 亚松弛
    void relax(Scalar alpha);

    // 是否有效（初始化与否）
    bool isValid() const;

private:
    // 私有求解辅助函数
    void ParallelSolve();
    void SerialSolve();


private:
    SparseMatrix<Tp>& equation_;  // 方程矩阵（稀疏矩阵，压不压缩都行）
    std::vector<Tp> x0_;        // 上一步的解
    std::vector<Tp> x_;         // 本步的解
    // Scalar relaxationFactor_{ 1.0 };// 松弛因子
    Method method_;             // 求解方法
    int maxIterationNum_;       // 最大迭代次数
    Scalar tolerance_;           // 残差
    Field<Tp>* filed_{ nullptr };           // 待求解的场
    // 每多少步输出一次
    int outputInterval_;
    bool isParallel_{ false };  // 是否并行
    bool isInitialized_{ false };// 是否赋初始值
};






template<typename Tp>
inline Solver<Tp>::Solver(SparseMatrix<Tp>& matrix, Method method, int maxmaxIterationNum)
    : equation_(matrix)
    , method_(method)
    , maxIterationNum_(maxmaxIterationNum)
    , tolerance_(DEFAULT_EPSILON)
    , outputInterval_(DEFAULT_OUTPUT_INTERVAL)
{
    // 检查矩阵是否有效
    if (!equation_.isValid())
    {
        std::cerr << "Solver<Tp>::Solver(const SparseMatrix<Tp>& matrix) Error: matrix is not valid" << std::endl;
        throw std::invalid_argument("matrix is not valid");
    }

    // 矩阵有效，判断是否压缩，没压缩，则进行压缩
    if (!equation_.isCompressed())
    {
        equation_.compress();
    }
}

template<typename Tp>
inline void Solver<Tp>::init(const std::vector<Tp>& x0)
{
    // 检查长度是否合法
    if (x0.size() != equation_.size())
    {
        std::cerr << "Solver<Tp>::init(const std::vector<Tp>& x0) Error: length of x0 is not equal to matrix size, length of x0i is illegal" << std::endl;
        throw std::invalid_argument("length of x0 is not equal to matrix size, length of x0i is illegal");
    }
    // 赋初值
    x0_ = x0;
    x_.resize(equation_.size());

    // 获取场，如果没有就是空指针
    filed_ = &equation_.getField();

    isInitialized_ = true;
}

template<typename Tp>
inline void Solver<Tp>::init(Tp value0)
{
    // 不能重复初始化
    if (isInitialized_)
    {
        std::cerr << "Solver<Tp>::init(Tp value0) Error: cannot initialize again" << std::endl;
        throw std::runtime_error("cannot initialize again");
    }

    x0_.resize(equation_.size(), value0);
    x_.resize(equation_.size());
    isInitialized_ = true;
}

template<typename Tp>
inline void Solver<Tp>::init(const Field<Tp>& initField)
{
    if (isInitialized_)
    {
        std::cerr << "Solver<Tp>::init(const Field<Tp>& initField) Error: cannot initialize again" << std::endl;
        throw std::runtime_error("cannot initialize again");
    }

    // 判断场是否有效
    if (!initField.isValid())
    {
        std::cerr << "Solver<Tp>::init(const Field<Tp>& initField) Error: initField is not valid" << std::endl;
        throw std::invalid_argument("initField is not valid");
    }

    Mesh* equationMesh = equation_.getMesh();
    Mesh* initFieldMesh = initField.getMesh();

    // 判断是否为同一个场
    // 如果稀疏矩阵不是由网格构造，则可忽略
    if (equationMesh != nullptr &&
        initFieldMesh != equationMesh)
    {
        std::cerr << "The mesh of the field is not equal to the mesh of the matrix" << std::endl;
        throw std::invalid_argument("The mesh of the field is not equal to the mesh of the matrix");
    }

    // 判断长度是否合理
    if (equation_.size() != initFieldMesh->getCellNumber())
    {
        std::cerr << "Solver<Tp>::init(const Field<Tp>& initField) Error: The size of equation matrix is not equal to the number of cells in the mesh" << std::endl;
        throw std::invalid_argument("The size of equation matrix is not equal to the number of cells in the mesh");
    }

    // 从field中的cellField_0_获取初始值
    x0_ = initField.getCellField_0().getData();
    x_.resize(equation_.size());
    isInitialized_ = true;
}

template<typename Tp>
inline Scalar Solver<Tp>::Error()
{
    // const std::vector<Scalar>& values = equation_.getValues();
    // const std::vector<ULL>& rowPointer = equation_.getRowPointer();
    // const std::vector<ULL>& columnIndex = equation_.getColIndexs();
    // const std::vector<Tp>& b = equation_.getB();

    ULL size = equation_.size();
    const Field<Tp>& phi = equation_.getField();

    // 获取本步解和上一部解
    const std::vector<Tp>& x = phi.getCellField().getData();    // 本步解
    const std::vector<Tp>& x0 = phi.getCellField_0().getData(); // 上一部解

    // 采用相对误差方式
    if constexpr (std::is_same_v<Tp, Scalar>)   // 标量场
    {
        Scalar temp1{};     // 分子
        Scalar temp2{};     // 分母

        for (ULL i = 0; i < size; ++i)
        {
            temp1 += (x[i] - x0[i]) * (x[i] - x0[i]);
            temp2 += x[i] * x[i];
        }

        if (temp1 == 0) // 第一次直接返回最大值
        {
            return std::numeric_limits<Scalar>::max();
        }

        return std::sqrt(temp1 / (temp2 + 1e-30));
    }
    else if constexpr (std::is_same_v<Tp, Vector<Scalar>>)  // 矢量场用相减后的模长平方
    {
        Scalar temp1{};
        Scalar temp2{};

        for (ULL i = 0; i < size; ++i)
        {
            temp1 += (x[i] - x0[i]).magnitudeSquared();
            temp2 += x[i].magnitudeSquared();
        }

        if (temp1 == 0) // 第一次直接返回最大值
        {
            return std::numeric_limits<Scalar>::max();
        }

        return std::sqrt(temp1 / (temp2 + 1e-30));
    }




    // if constexpr (std::is_same_v<Tp, Scalar>)
    // {
    //     for (ULL i = 0; i < size; ++i)
    //     {
    //         Scalar residual{};
    //         for (ULL colId = rowPointer[i]; colId < rowPointer[i + 1]; ++colId)
    //         {
    //             residual += values[colId] * x0_[columnIndex[colId]];
    //         }
    //         residual = std::abs(b[i] - residual);
    //         maxResidual = std::max(maxResidual, residual);
    //     }
    // }
    // else if constexpr (std::is_same_v<Tp, Vector<Scalar>>)
    // {
    //     for (ULL i = 0; i < size; ++i)
    //     {
    //         Vector<Scalar> residual{};
    //         for (ULL colId = rowPointer[i]; colId < rowPointer[i + 1]; ++colId)
    //         {
    //             residual += values[colId] * x0_[columnIndex[colId]];
    //         }
    //         residual = residual - b[i];
    //         maxResidual = std::max(maxResidual, residual.magnitude());
    //     }
    // }
    // else
    // {
    //     std::cerr << "Solver<Tp>::SerialSolve() Error: The type of Tp is not supported" << std::endl;
    //     throw std::runtime_error("The type of Tp is not supported");
    // }
}



template<typename Tp>
inline void Solver<Tp>::setParallel()
{
    isParallel_ = true;
}

template<typename Tp>
inline void Solver<Tp>::solve()
{
    if (!isInitialized_)
    {
        std::cerr << "Solver<Tp>::solve() Error: cannot solve without initialization" << std::endl;
        throw std::runtime_error("cannot solve without initialization");
    }

    if (isParallel_)
    {
        ParallelSolve();
    }
    else
    {
        SerialSolve();
    }


}

template<typename Tp>
inline void Solver<Tp>::setTolerance(Scalar tolerance)
{
    if (!isInitialized_)
    {
        std::cerr << "Solver<Tp>::setTolerance() Error: cannot set tolerance without initialization" << std::endl;
        throw std::runtime_error("cannot set tolerance without initialization");
    }
    tolerance_ = tolerance;
}

template<typename Tp>
inline void Solver<Tp>::relax(Scalar alpha)
{
    if (!isInitialized_)
    {
        std::cerr << "Solver<Tp>::relax() Error: cannot relax without initialization" << std::endl;
        throw std::runtime_error("cannot relax without initialization");
    }

    // 判断alpha是否在0-1之间
    if (alpha <= 0 || alpha >= 1)
    {
        std::cerr << "Solver<Tp>::relax() Error: alpha must be in (0, 1)" << std::endl;
        throw std::runtime_error("alpha must be in (0, 1)");
    }



    // 将方程对角线元素除以alpha
    ULL size = equation_.size();
    std::vector<Scalar> diag(size);
    for (ULL i = 0; i < size; ++i)
    {
        // 对角线记录的是已经除以alpha后的值
        diag[i] = (equation_(i, i) /= alpha);
    }

    // 右端项
    for (ULL i = 0; i < size; ++i)
    {
        equation_.addB(i, (1 - alpha) * diag[i] * x0_[i]);
    }
}

template<typename Tp>
inline bool Solver<Tp>::isValid() const
{
    return isInitialized_;
}

template<typename Tp>
inline void Solver<Tp>::ParallelSolve()
{
    if (method_ != Solver::Method::Jacobi)
    {
        std::cerr << "Solver<Tp>::ParallelSolve() Error: only Jacobi method is supported"
            << std::endl;
        throw std::invalid_argument("only Jacobi method is supported");
    }

    const ULL size = equation_.size();

    const std::vector<ULL>& rowPointer = equation_.getRowPointer();
    const std::vector<ULL>& columnIndex = equation_.getColIndexs();
    const std::vector<Scalar>& values = equation_.getValues();
    const std::vector<Tp>& b = equation_.getB();

    // 保存求解前的初始解。
    // 这样最后写回 Field 时，cellField_0 可以表示本次线性求解前的旧场。
    const std::vector<Tp> oldSolution = x0_;

    // 提前提取对角线，避免每次迭代都调用 equation_(i, i) 查找 CSR。
    std::vector<Scalar> diag(size, Scalar{});

    for (ULL row = 0; row < size; ++row)
    {
        if (rowPointer[row] == rowPointer[row + 1])
        {
            std::cerr << "Solver<Tp>::ParallelSolve() Error: row "
                << row << " is an empty row" << std::endl;
            throw std::runtime_error("matrix has an empty row");
        }

        bool hasDiagonal = false;

        for (ULL index = rowPointer[row]; index < rowPointer[row + 1]; ++index)
        {
            if (columnIndex[index] == row)
            {
                diag[row] = values[index];
                hasDiagonal = true;
                break;
            }
        }

        if (!hasDiagonal || std::abs(diag[row]) < static_cast<Scalar>(1.0e-30))
        {
            std::cerr << "Solver<Tp>::ParallelSolve() Error: invalid diagonal at row "
                << row << std::endl;
            throw std::runtime_error("invalid diagonal element");
        }
    }

    // 线程数不是越多越好。
    // Jacobi + CSR SpMV 通常受内存带宽限制，矩阵太小时不值得开太多线程。
    const unsigned hardwareThreadNum = std::thread::hardware_concurrency();
    const unsigned maxThreadNum = hardwareThreadNum > 1U ? hardwareThreadNum - 1U : 1U;

    constexpr ULL MIN_ROWS_PER_TASK = 1024;

    const ULL usefulThreadNum = std::max<ULL>(
        1ULL,
        (size + MIN_ROWS_PER_TASK - 1ULL) / MIN_ROWS_PER_TASK
    );

    const unsigned threadNum = static_cast<unsigned>(
        std::min<ULL>(
            static_cast<ULL>(maxThreadNum),
            usefulThreadNum
        )
        );

    if (threadNum <= 1U)
    {
        SerialSolve();
        return;
    }

    struct Range
    {
        ULL begin{};
        ULL end{};
    };

    std::vector<Range> ranges(threadNum);

    const ULL baseLineNum = size / threadNum;
    const ULL remainLineNum = size % threadNum;

    ULL begin = 0;

    for (unsigned tid = 0; tid < threadNum; ++tid)
    {
        const ULL lineNum = baseLineNum + (tid < remainLineNum ? 1ULL : 0ULL);

        ranges[tid] = Range{
            begin,
            begin + lineNum
        };

        begin += lineNum;
    }

    ThreadPool& pool = ThreadPool::getInstance();

    // 如果线程池已经启动，新的 ThreadPool::start() 会直接 return，不会重复创建线程。
    pool.start(threadNum);

    std::vector<std::future<void>> updateFutures(threadNum);
    std::vector<std::future<Scalar>> residualFutures(threadNum);

    // 第一阶段：Jacobi 更新。
    //
    // 每个线程只写自己的连续区间：
    //
    //     x_[rowBegin, rowEnd)
    //
    // 因此不会有多个线程写同一个 x_[i]。
    auto updateRange = [&](ULL rowBegin, ULL rowEnd) {
        for (ULL row = rowBegin; row < rowEnd; ++row)
        {
            Tp rhs = b[row];

            // Jacobi 迭代：
            //
            //     x_i^{k+1}
            //       =
            //     ( b_i - sum_{j != i} A_ij x_j^k ) / A_ii
            //
            // 这里只读旧解 x0_，只写新解 x_。
            for (ULL index = rowPointer[row]; index < rowPointer[row + 1]; ++index)
            {
                const ULL col = columnIndex[index];

                if (col == row)
                {
                    continue;
                }

                rhs -= values[index] * x0_[col];
            }

            x_[row] = rhs / diag[row];
        }
        };

    // 第二阶段：残差计算。
    //
    // 每个线程返回自己的局部最大残差，主线程最后做 max 归约。
    // 不使用 mutex，也不使用 atomic。
    auto residualRange = [&](ULL rowBegin, ULL rowEnd) -> Scalar {
        Scalar localMaxResidual{};

        if constexpr (std::is_same_v<Tp, Scalar>)
        {
            for (ULL row = rowBegin; row < rowEnd; ++row)
            {
                Scalar ax{};

                for (ULL index = rowPointer[row]; index < rowPointer[row + 1]; ++index)
                {
                    ax += values[index] * x_[columnIndex[index]];
                }

                const Scalar residual = std::abs(b[row] - ax);
                localMaxResidual = std::max(localMaxResidual, residual);
            }
        }
        else if constexpr (std::is_same_v<Tp, Vector<Scalar>>)
        {
            for (ULL row = rowBegin; row < rowEnd; ++row)
            {
                Vector<Scalar> ax{};

                for (ULL index = rowPointer[row]; index < rowPointer[row + 1]; ++index)
                {
                    ax += values[index] * x_[columnIndex[index]];
                }

                const Scalar residual = (b[row] - ax).magnitude();
                localMaxResidual = std::max(localMaxResidual, residual);
            }
        }
        else
        {
            std::cerr << "Solver<Tp>::ParallelSolve() Error: unsupported Tp type"
                << std::endl;
            throw std::runtime_error("unsupported Tp type");
        }

        return localMaxResidual;
        };

    for (int it = 0; it < maxIterationNum_; ++it)
    {
        // 并行更新 x_。
        for (unsigned tid = 0; tid < threadNum; ++tid)
        {
            updateFutures[tid] = pool.submitTask(
                updateRange,
                ranges[tid].begin,
                ranges[tid].end
            );
        }

        for (unsigned tid = 0; tid < threadNum; ++tid)
        {
            updateFutures[tid].get();
        }

        // 并行计算残差。
        // 必须等所有 x_ 更新完成后才能计算 residual。
        for (unsigned tid = 0; tid < threadNum; ++tid)
        {
            residualFutures[tid] = pool.submitTask(
                residualRange,
                ranges[tid].begin,
                ranges[tid].end
            );
        }

        Scalar maxResidual{};

        for (unsigned tid = 0; tid < threadNum; ++tid)
        {
            maxResidual = std::max(maxResidual, residualFutures[tid].get());
        }


        if (maxResidual < tolerance_ || it == maxIterationNum_ - 1)
        {
            // 写回 Field。
            //
            // cellField    ：本次线性求解后的新解
            // cellField_0  ：本次线性求解前的旧解
            if (filed_ != nullptr)
            {
                filed_->getCellField().getData() = x_;
                filed_->getCellField_0().getData() = oldSolution;
            }

            return;
        }

        // Jacobi 下一轮只需要交换新旧解。
        // 不要写 x0_ = x_，那会整数组拷贝。
        x0_.swap(x_);
    }
}

template<typename Tp>
inline void Solver<Tp>::SerialSolve()
{
    // 获取必要参数
    ULL size = equation_.size();    // size
    const std::vector<ULL>& rowPointer = equation_.getRowPointer(); // 稀疏矩阵行首索引
    const std::vector<ULL>& columnIndex = equation_.getColIndexs(); // 与values对应的列索引
    const std::vector<Scalar>& values = equation_.getValues();   // 稀疏非0元素
    const std::vector<Tp>& b = equation_.getB();


    // 先确认不存在0元素行
    for (ULL i = 0; i < size; ++i)
    {
        if (rowPointer[i] == rowPointer[i + 1])
        {
            std::cerr << "Solver<Tp>::SerialSolve() Error: The line " << i << " of the matrix is a zero element row" << std::endl;
            throw std::runtime_error("matrix has a zero element row");
        }
    }

    // 获取矩阵对角元素
    std::vector<Scalar> diag(size);
    for (ULL i = 0; i < size; ++i)
    {
        diag[i] = equation_(i, i);
    }

    if (method_ == Solver::Method::Jacobi)
    {
        for (int it = 0; it < maxIterationNum_; ++it)
        {
            for (ULL i = 0; i < size; ++i)
            {
                Tp sum{};

                // 只遍历当前行的非0元素
                for (ULL colId = rowPointer[i]; colId < rowPointer[i + 1]; ++colId)
                {
                    sum += values[colId] * x0_[columnIndex[colId]];
                }
                sum -= x0_[i] * diag[i];
                // 计算当前步解
                x_[i] = (b[i] - sum) / diag[i];
            }
            // x0_ = x_;   // 更新x0_

            // 计算残差，采用最大残差
            Scalar maxResidual = 0;
            if constexpr (std::is_same_v<Tp, Scalar>)
            {
                for (ULL i = 0; i < size; ++i)
                {
                    Scalar residual{};
                    for (ULL colId = rowPointer[i]; colId < rowPointer[i + 1]; ++colId)
                    {
                        residual += values[colId] * x_[columnIndex[colId]];
                    }
                    residual = std::abs(b[i] - residual);
                    maxResidual = std::max(maxResidual, residual);
                }
            }
            else if constexpr (std::is_same_v<Tp, Vector<Scalar>>)
            {
                for (ULL i = 0; i < size; ++i)
                {
                    Vector<Scalar> residual{};
                    for (ULL colId = rowPointer[i]; colId < rowPointer[i + 1]; ++colId)
                    {
                        residual += values[colId] * x_[columnIndex[colId]];
                    }
                    residual = residual - b[i];
                    maxResidual = std::max(maxResidual, residual.magnitude());
                }
            }
            else
            {
                std::cerr << "Solver<Tp>::SerialSolve() Error: The type of Tp is not supported" << std::endl;
                throw std::runtime_error("The type of Tp is not supported");
            }

            // // 每N步输出一次
            // if ((it + 1) % outputInterval_ == 0 ||  // 保证后续输出的是整数次
            //     it % outputInterval_ == 0)  // 第0次用
            // {
            //     if (filed_)
            //     {
            //         filed_->getCellField().getData() = x_;
            //         filed_->getCellField_0().getData() = x_;
            //     }
            //     std::cout << "Iteration " << it + 1 << ": Max Residual = " << maxResidual << std::endl;
            // }



            // 测试用
            // for (int i = 0; i < equation_.size(); ++i)
            // {

            //     std::cout << x_[i] << " ";
            // }
            // std::cout << std::endl;

            // 残差满足要求或到最大迭代次数则返回，本次求解结束
            if (maxResidual < tolerance_ ||
                it == maxIterationNum_ - 1)
            {
                if (filed_)
                {
                    // filed_->getCellField().getData() = x_;

                    // 新值给cellField_0的，再与cellField交换
                    filed_->getCellField_0().getData() = std::move(x_);
                    filed_->getCellField_0().getData().
                        swap(filed_->getCellField().getData());
                }
                // std::cout << "Iteration " << it + 1 << ": Max Residual = " << maxResidual << " The solution has converged" << std::endl;

                return;
            }

            // 更新x0_
            x0_.swap(x_);
        }
    }
    else if (method_ == Solver::Method::GaussSeidel)
    {

    }
    else if (method_ == Solver::Method::AMG)
    {

    }
}






#endif // SOLVER_H_