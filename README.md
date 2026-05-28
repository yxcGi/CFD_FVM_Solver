# 摘要

**CFD_FVM_Solver** 是一个基于 C++ 开发的轻量级计算流体力学（CFD）求解器。该项目采用**有限体积法（Finite Volume Method, FVM）**作为核心离散算法，旨在解决不可压缩流体的对流扩散问题。

该求解器采用面向对象的设计思想（OOP），实现了非结构网格的读取与处理、标量/矢量场的管理、稀疏矩阵存储与求解以及并行计算支持。项目结构类似于 OpenFOAM，但在底层实现上更加轻量且易于理解，适合用于 CFD 算法研究与教学。

---
# ⭐️Recent Updates
1. 添加了SIMPLE算法。
2. 修复了二阶迎风(SUD)延迟修正系数的错误。

# 主要特性

## 1. 网格与几何 (Mesh & Geometry)
*   **非结构网格支持**：支持任意多面体网格（包括三角形、四边形、四面体、六面体等）。
*   **OpenFOAM 格式兼容**：直接支持读取 OpenFOAM 的 `polyMesh` 网格格式（`points`, `faces`, `owner`, `neighbour`, `boundary`）。
*   **几何计算**：自动计算单元体积、面面积、面法向量、单元中心及面中心。
*   **边界类型**：支持 `patch` (进出口), `wall` (壁面), `symmetry` (对称), `empty` (二维模拟), `cyclic` (周期性) 等多种边界类型。

## 2. 场与数据结构 (Fields & Math)
*   **场操作**：实现了 `CellField` (体心场) 和 `FaceField` (面心场)，支持标量（Scalar）与矢量（Vector）场。
*   **边界条件**：支持混合边界条件配置 ($a\phi + b\frac{\partial \phi}{\partial n} = c$)，可轻松实现第一类（Dirichlet）、第二类（Neumann）及 Robin 边界条件。
*   **数学库**：内置高性能的 `Vector` 和 `Tensor` 模板类，支持各种张量运算。
*   **稀疏矩阵**：采用 **CSR (Compressed Sparse Row)** 格式存储稀疏矩阵，优化内存占用与访问速度。

## 3. 离散格式 (Discretization)
*   **对流项 (Convection)**：
    *   一阶迎风格式 (FUD)
    *   二阶迎风格式 (SUD，带延迟修正)
    *   中心差分格式 (CD)
*   **扩散项 (Laplacian)**：包含非正交修正（Non-orthogonal correction），适用于扭曲网格。
*   **梯度计算**：基于高斯-格林公式（Gauss-Green）的梯度重构。
*   **插值**：支持线性插值、算术平均、调和平均等多种插值方案。

## 4. 线性方程组求解 (Linear Solver)
*   **求解器**：内置 Jacobi 迭代求解器。
*   **并行计算**：基于自定义 **线程池 (Thread Pool)** 实现的多线程共享内存并行计算，加速矩阵求解过程。
*   **残差控制**：支持自定义收敛容差与最大迭代次数。

## 5. 后处理 (Post-processing)
*   **Tecplot 输出**：支持将计算结果输出为 Tecplot (`.dat`) 格式，兼容 2D 和 3D 数据的可视化。

---

# 目录结构

```text
src
├── CMakeLists.txt                      # 顶层 CMake 构建脚本
├── core                                # 核心库目录（包含求解器的底层实现）
│   ├── CMakeLists.txt                  # 核心库的 CMake 构建脚本
│   ├── FVM                             # 有限体积法离散算子（Discretization）
│   │   ├── Div.hpp                     # 对流项（散度）离散函数
│   │   ├── DivType.h                   # 对流项离散格式枚举（如一阶迎风、二阶迎风、中心差分等）
│   │   └── Laplacian.hpp               # 扩散项（拉普拉斯项）离散函数
│   ├── Fields                          # 物理场数据结构与管理
│   │   ├── BaseField.hpp               # 场类的基类（定义通用接口）
│   │   ├── BoundaryCondition.hpp       # 边界条件类（处理第一、二、三类边界条件）
│   │   ├── CellField.hpp               # 体心场（存储在单元中心的数据）
│   │   ├── FaceField.hpp               # 面心场（存储在面中心的数据）
│   │   ├── Field.hpp                   # 场的高级封装（包含新旧时间步数据、梯度等）
│   │   ├── FieldOperators              # 场运算算子
│   │   │   ├── Divergence.hpp          # 场的散度计算
│   │   │   ├── Gradient.hpp            # 场的梯度计算（如高斯-格林公式）
│   │   │   └── Laplacian.hpp           # 场的拉普拉斯计算
│   │   └── FieldType.hpp               # 场类型枚举
│   ├── Geometry                        # 网格几何拓扑处理
│   │   ├── BoundaryPatch.cpp           # 边界补丁实现
│   │   ├── BoundaryPatch.h             # 边界补丁定义（定义边界类型、起止面索引）
│   │   ├── Cell.cpp                    # 网格单元实现（体积、中心计算）
│   │   ├── Cell.h                      # 网格单元定义
│   │   ├── Face.cpp                    # 网格面实现（面积、法向量计算）
│   │   ├── Face.h                      # 网格面定义
│   │   ├── Mesh.cpp                    # 网格整体管理实现（读取 OpenFOAM 格式网格）
│   │   └── Mesh.h                      # 网格类定义
│   ├── Math                            # 数学基础库与线性代数求解器
│   │   ├── Solver.hpp                  # 线性方程组求解器（如 Jacobi 迭代）
│   │   ├── SparseMatrix.hpp            # 稀疏矩阵类（CSR 格式存储）
│   │   ├── Tensor.hpp                  # 二阶张量类及运算
│   │   ├── Vector.hpp                  # 三维向量类及运算
│   │   └── interpolation.hpp           # 插值函数（从体心到面心的插值方案）
│   └── thread_pool                     # 并行计算工具
│       └── threadpool.hpp              # 线程池实现（用于矩阵求解并行化）
└── main.cpp                            # 程序主入口（包含测试案例与求解流程）
```

---

# 编译与构建

本项目使用 **CMake** 进行构建管理。请确保您的环境中已安装 C++17 编译器（如 GCC, Clang, MSVC）和 CMake，推荐使用Linux 或 macOS 系统。

1.  **克隆或下载代码**
    ```bash
    git clone <repository_url>
    cd CFD_FVM_Solver
    ```

2.  **创建构建目录并编译**
    ```bash
    cmake -B build
    cmake --build build
    ```

3.  **运行**
    ```bash
    ./bin/CFD_FVM_Solver
    ```

---

# 使用示例

`main.cpp` 中有若干个对流扩散的案例可供参考。

在 `main.cpp` 中可以配置具体的物理问题。以下是一个简单的热传导（扩散）或对流扩散问题的配置流程：

### 1. 读取网格
```cpp
// 指定 OpenFOAM polyMesh 文件夹路径
Mesh mesh("path/to/your/cavity/constant/polyMesh");
```

### 2. 创建物理场
```cpp
// 创建温度场 T
Field<Scalar> T("T", &mesh);
T.setValue(300); // 初始化为 300K

// 设置边界条件
// setBoundaryCondition(边界名, a, b, c) -> a*T + b*dT/dn = c
// 例如：Dirichlet (T=500): a=1, b=0, c=500
//       Neumann (dT/dn=0): a=0, b=1, c=0
T.setBoundaryCondition("movingWall", 1, 0, 500); 
T.setBoundaryCondition("fixedWalls", 1, 0, 300);
T.cellToFace(); // 更新边界值
```

### 3. 定义方程与求解
```cpp
// 定义扩散系数
FaceField<Scalar> gamma("gamma", &mesh);
gamma.setValue(0.01);

// 构建稀疏矩阵
SparseMatrix<Scalar> A(&mesh);

// 时间步/迭代循环
for (int i = 0; i < 1000; i++)
{
    // 离散化：离散扩散项
    fvm::Laplacian(A, gamma, T);

    // 求解线性方程组
    Solver<Scalar> solver(A, Solver<Scalar>::Method::Jacobi, 1000);
    solver.setParallel(); // 开启并行
    solver.init(T.getCellField_0().getData());
    solver.solve();

    // 更新场数据
    // ... (数据回写逻辑)
    
    // 清空矩阵用于下一次迭代
    A.clear();
}
```

### 4. 输出结果
```cpp
// 输出为 Tecplot 格式
T.writeToFile("Result_T.dat");
```

---

# 开发计划

*   [ ] 完善动量方程求解（SIMPLE/PISO 算法实现）。
*   [ ] 增加 AMG (代数多重网格) 求解器。
*   [ ] 支持瞬态项离散 (时间步进)。
*   [ ] 优化大规模网格下的内存管理。

---
# 联系作者
邮箱：<EMAIL>yuanxch3@mail2.sysu.edu.cn
