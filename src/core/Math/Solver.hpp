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
    // 不能重复初始化
    if (isInitialized_)
    {
        std::cerr << "Solver<Tp>::init(const std::vector<Tp>& x0) Error: cannot initialize again" << std::endl;
        throw std::runtime_error("cannot initialize again");
    }

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
    // 根据电脑cpu核数，创建线程池个数为 cpu核数-1，主线程用于集中处理
    ThreadPool* pool = &ThreadPool::getInstance();
    int threadNum = std::thread::hardware_concurrency() - 1;
    pool->start(threadNum);

    // 给每个线程分配计算行的范围
    ULL baseLineNum = equation_.size() / threadNum;
    ULL remainLineNum = equation_.size() % threadNum;
    std::vector<ULL> lineRange(threadNum + 1);  // 存储每个线程的行范围
    for (int i = 0; i < threadNum; ++i)     // 给范围赋值
    {
        if (i < remainLineNum)
        {
            lineRange[i] = (baseLineNum + 1) * i;
            lineRange[i + 1] = (baseLineNum + 1) * (i + 1);
        }
        else
        {
            lineRange[i] = baseLineNum * i + remainLineNum;
            lineRange[i + 1] = baseLineNum * (i + 1) + remainLineNum;
        }
    }

    ULL size = equation_.size();    // size
    const std::vector<ULL>& rowPointer = equation_.getRowPointer();
    const std::vector<Scalar>& values = equation_.getValues();
    const std::vector<ULL>& columnIndex = equation_.getColIndexs();
    const std::vector<Tp>& b = equation_.getB();
    std::vector<Scalar> diag(size);
    for (ULL i = 0; i < size; ++i)
    {
        diag[i] = equation_(i, i);
    }


    if (method_ == Solver::Method::Jacobi)
    {
        // 创建线程函数
        std::function<void(ULL, ULL)> paraSolve = [&](ULL start, ULL end) {
            for (ULL i = start; i < end; ++i)
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
            };

        // std::function<Scalar(ULL, ULL)> paraSolveFunc = std::bind(
        //     static_cast<Scalar(Solver<Tp>::*)(ULL, ULL)>(&Solver<Tp>::ParallelSolve),
        //     this,
        //     std::placeholders::_1,
        //     std::placeholders::_2
        // );

        // 存储线程的返回值数组
        std::vector<std::future<void>> results(threadNum);    // 记录每个线程返回的最大残差
        // auto start = std::chrono::high_resolution_clock::now();
        for (int it = 0; it < maxIterationNum_; ++it)
        {
            for (int i = 0; i < threadNum; ++i) // 给每个线程分配任务
            {
                results[i] = pool->submitTask(
                    paraSolve,
                    lineRange[i],
                    lineRange[i + 1]
                );
            }
            for (int i = 0; i < threadNum; ++i) // 等待线程完成
            {
                results[i].get();
            }

            x0_ = x_;   // 更新x0_

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

            // 每N步输出一次
            if ((it + 1) % outputInterval_ == 0 ||  // 保证后续输出的是整数次
                it % outputInterval_ == 0)  // 第0次用
            {
                if (filed_)
                {
                    filed_->getCellField().getData() = x_;
                    filed_->getCellField_0().getData() = x0_;
                }
                std::cout << "Iteration " << it + 1 << ": Max Residual = " << maxResidual << std::endl;
            }

            // 最大残差小于tolerance_，则结束迭代
            if (maxResidual < tolerance_)
            {
                if (filed_)
                {
                    filed_->getCellField().getData() = x_;
                    filed_->getCellField_0().getData() = x0_;
                }
                // std::cout << "Iteration " << it + 1 << ": Max Residual = " << maxResidual << "The solution has converged" << std::endl;

                // auto end = std::chrono::high_resolution_clock::now();
                // std::cout << "计算耗时：" << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;

                return;
            }
        }
    }
    else
    {
        std::cerr << "Solver<Tp>::SerialSolve() Error: The method is not supported" << std::endl;
        throw std::invalid_argument("The method is not supported");
    }
}



template<typename Tp>
inline void Solver<Tp>::SerialSolve()
{
    // 获取必要参数
    ULL size = equation_.size();    // size
    const std::vector<ULL>& rowPointer = equation_.getRowPointer(); // 稀疏矩阵行首索引
    const std::vector<ULL>& columnIndex = equation_.getColIndexs(); // 与values对应的列索引
    const std::vector<Scalar> values = equation_.getValues();   // 稀疏非0元素
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