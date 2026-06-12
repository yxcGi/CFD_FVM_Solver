# CFD_FVM_Solver

**CFD_FVM_Solver** 是一个基于 **C++20** 开发的轻量级有限体积法（Finite Volume Method, FVM）CFD 求解器。项目目标是从底层实现网格读取、几何计算、物理场管理、有限体积离散、稀疏矩阵装配、线性方程组求解以及 SIMPLE 不可压缩流求解流程。

该项目的设计思路接近 OpenFOAM 的基本抽象方式：

- `Mesh`：网格与几何拓扑；
- `Field`：物理场；
- `fvm::Laplacian / fvm::Div / fvm::Source`：有限体积离散算子；
- `SparseMatrix`：稀疏矩阵装配；
- `Solver`：线性方程组迭代求解；
- `algorithm::simple::SIMPLE`：不可压缩稳态流 SIMPLE 算法；
- `example/`：独立算例管理目录。

项目主要用于 **CFD 算法学习、有限体积法开发、非结构网格离散验证、SIMPLE 算法实现验证**。

---

# 算例展示

以下算例覆盖二维方腔驱动、二维台阶流与三维方腔驱动等典型不可压缩流验证场景。

## 二维方腔驱动

方腔驱动算例同时展示了三角形非结构化网格与正方形结构化网格下的流场结果。

<table>
  <tr>
    <td align="center">
      <img src="docs/example_image/cavity2d_triangle_unstructured_mesh.png" width="280"><br>
      三角形非结构化网格
    </td>
    <td align="center">
      <img src="docs/example_image/cavity2d_square_structured_mesh_20x20.png" width="280"><br>
      正方形结构化网格 20×20
    </td>
  </tr>
</table>

| Re | 正方形网格 | 三角形网格 |
| --- | --- | --- |
| 400 | <img src="docs/example_image/cavity2d_square_re400.png" width="220"> | <img src="docs/example_image/cavity2d_triangle_re400.png" width="220"> |
| 1000 | <img src="docs/example_image/cavity2d_square_re1000.png" width="220"> | <img src="docs/example_image/cavity2d_triangle_re1000.png" width="220"> |
| 2000 | <img src="docs/example_image/cavity2d_square_re2000.png" width="220"> | <img src="docs/example_image/cavity2d_triangle_re2000.png" width="220"> |
| 5000 | <img src="docs/example_image/cavity2d_square_re5000.png" width="220"> | <img src="docs/example_image/cavity2d_triangle_re5000.png" width="220"> |

## 二维台阶流

二维台阶流算例网格数为 12225，结果展示了不同 Reynolds 数下的流线与速度云图。

<p align="center">
  <img src="docs/example_image/backward_step_mesh_12225.png" width="720"><br>
  台阶流计算网格
</p>

| Re | 流线 | 速度云图 |
| --- | --- | --- |
| 50 | <img src="docs/example_image/backward_step_streamline_re50.png" width="320"> | <img src="docs/example_image/backward_step_velocity_re50.png" width="320"> |
| 100 | <img src="docs/example_image/backward_step_streamline_re100.png" width="320"> | <img src="docs/example_image/backward_step_velocity_re100.png" width="320"> |
| 150 | <img src="docs/example_image/backward_step_streamline_re150.png" width="320"> | <img src="docs/example_image/backward_step_velocity_re150.png" width="320"> |
| 200 | <img src="docs/example_image/backward_step_streamline_re200.png" width="320"> | <img src="docs/example_image/backward_step_velocity_re200.png" width="320"> |

## 三维方腔驱动

三维方腔驱动算例网格数为 10250，边长 `L = 1`，下表展示了不同截面和 Reynolds 数下的速度场结构。

<table>
  <tr>
    <td align="center">
      <img src="docs/example_image/cavity3d_velocity_direction_overview.png" width="420"><br>
      三维方腔速度方向示意
    </td>
    <td align="center">
      <img src="docs/example_image/cavity3d_velocity_direction_axes.png" width="320"><br>
      坐标轴与速度方向
    </td>
  </tr>
</table>

| 截面 | Re = 400 | Re = 1000 | Re = 1600 |
| --- | --- | --- | --- |
| z = 0.5 | <img src="docs/example_image/cavity3d_z05_re400.png" width="180"> | <img src="docs/example_image/cavity3d_z05_re1000.png" width="180"> | <img src="docs/example_image/cavity3d_z05_re1600.png" width="180"> |
| x = 0.5 | <img src="docs/example_image/cavity3d_x05_re400.png" width="180"> | <img src="docs/example_image/cavity3d_x05_re1000.png" width="180"> | <img src="docs/example_image/cavity3d_x05_re1600.png" width="180"> |
| y = 0.5 | <img src="docs/example_image/cavity3d_y05_re400.png" width="180"> | <img src="docs/example_image/cavity3d_y05_re1000.png" width="180"> | <img src="docs/example_image/cavity3d_y05_re1600.png" width="180"> |

---

## 1. 当前功能

### 1.1 网格与几何模块

当前支持读取 OpenFOAM 风格的 `polyMesh` 网格目录：

```text
constant/polyMesh/
├── points
├── faces
├── owner
├── neighbour
└── boundary
```

`Mesh` 会管理：

- 点坐标 `points`；
- 面拓扑 `faces`；
- owner / neighbour 单元关系；
- 内部面索引；
- 边界面索引；
- 单元列表；
- 边界 patch 信息；
- 网格维度；
- Tecplot 输出所需的网格形状信息。

支持的边界 patch 类型包括：

```cpp
PATCH
WALL
SYMMETRY
CYCLIC
WEDGE
EMPTY
PROCESSOR
CUSTOM
```

其中 `EMPTY` 主要用于二维算例中的空边界。

---

### 1.2 几何实体

项目中主要有三个几何实体：

```text
Mesh
├── Face
├── Cell
└── BoundaryPatch
```

#### Face

`Face` 负责保存和计算：

- 面的点索引；
- owner 单元索引；
- neighbour 单元索引；
- 面心；
- 面积；
- 单位法向量。

#### Cell

`Cell` 负责保存和计算：

- 单元包含的面索引；
- 单元包含的点索引；
- 单元中心；
- 单元体积。

#### BoundaryPatch

`BoundaryPatch` 负责记录：

- patch 名称；
- patch 类型；
- patch 起始面索引；
- patch 面数量。

---

### 1.3 场变量系统

项目中的物理场主要由以下类组成：

```text
BaseField<Tp>
├── CellField<Tp>
├── FaceField<Tp>
└── Field<Tp>
```

其中：

- `CellField<Tp>`：体心场；
- `FaceField<Tp>`：面心场；
- `Field<Tp>`：高级物理场封装，同时包含体心场、面心场、旧体心场和梯度场。

常见用法：

```cpp
Mesh mesh("tempFile/OpenFOAM_tutorials/cavity/constant/polyMesh");

Field<Scalar> p("p", &mesh);
Field<Vector<Scalar>> U("U", &mesh);

p.setValue(0.0);
U.setValue(Vector<Scalar>(0, 0, 0));
```

也支持通过坐标函数初始化：

```cpp
Field<Scalar> T("T", &mesh);

T.setValue(
    [](Scalar x, Scalar y, Scalar z) {
        return x + y + z;
    }
);
```

---

### 1.4 边界条件

边界条件采用统一的 Robin 形式：

```math
a \phi + b \frac{\partial \phi}{\partial n} = c
```

对应代码接口：

```cpp
field.setBoundaryCondition("patchName", a, b, c);
```

常见边界条件写法：

| 类型 | 数学形式 | 代码写法 |
|---|---|---|
| Dirichlet | `phi = c` | `a = 1, b = 0` |
| Neumann | `dphi/dn = c` | `a = 0, b = 1` |
| Robin | `a phi + b dphi/dn = c` | 直接给定 `a, b, c` |

例如：

```cpp
T.setBoundaryCondition("leftWalls", 1, 0, 100.0);
T.setBoundaryCondition("rightWalls", 0, 1, 0.0);
T.cellToFace();
```

---

### 1.5 数学基础类

项目实现了基础数学类型：

```text
Vector<Tp>
Tensor<Tp>
```

#### Vector

`Vector<Tp>` 支持：

- 向量加减；
- 标量乘除；
- 点乘；
- 叉乘；
- 模长；
- 模长平方；
- 单位向量；
- 两点距离；
- 流输出。

示例：

```cpp
Vector<Scalar> a(1.0, 2.0, 3.0);
Vector<Scalar> b(4.0, 5.0, 6.0);

Scalar dot = a & b;
Vector<Scalar> cross = a ^ b;
Scalar mag = a.magnitude();
```

#### Tensor

`Tensor<Tp>` 支持：

- 二阶张量构造；
- 2D / 3D 构造；
- 对角初始化；
- 张量加减；
- 张量乘法；
- 张量与向量乘法；
- 双点积；
- 标量乘除。

---

### 1.6 插值模块

插值模块位于：

```text
src/core/Math/interpolation.hpp
```

当前支持：

```cpp
interpolation::Scheme::LINEAR
interpolation::Scheme::ARITHMETIC_MEAN
interpolation::Scheme::HARMONIC_MEAN
```

标量、向量、张量均有对应插值接口。常见调用方式：

```cpp
Scalar phiF = util::interpolate(
    phiOwner,
    phiNeighbor,
    interpolation::Scheme::LINEAR,
    alpha
);
```

对于扩散系数等物性参数，可以使用调和平均：

```cpp
Scalar gammaF = util::interpolate(
    gammaOwner,
    gammaNeighbor,
    interpolation::Scheme::HARMONIC_MEAN,
    alpha
);
```

---

## 2. 有限体积离散模块

有限体积离散算子位于：

```text
src/core/FVM/
├── Div.hpp
├── DivType.h
├── Laplacian.hpp
└── Source.hpp
```

---

### 2.1 扩散项 `fvm::Laplacian`

扩散项接口：

```cpp
fvm::Laplacian(matrix, gamma, phi);
```

其中：

- `matrix`：待装配的稀疏矩阵；
- `gamma`：面心扩散系数；
- `phi`：待离散物理场。

当前扩散项实现包含：

- 内部面扩散通量装配；
- 非正交修正；
- 显式延迟修正；
- Robin 边界条件处理；
- 标量场和矢量场统一模板处理。

典型用法：

```cpp
Field<Scalar> T("T", &mesh);
FaceField<Scalar> gamma("gamma", &mesh);
SparseMatrix<Scalar> A(&mesh);

T.setValue(500.0);
gamma.setValue(20.0);

T.setBoundaryCondition("leftWalls", 1, 0, 100.0);
T.setBoundaryCondition("rightWalls", 0, 1, 0.0);
T.cellToFace();

fvm::Laplacian(A, gamma, T);
```

---

### 2.2 对流项 `fvm::Div`

对流项接口：

```cpp
fvm::Div(matrix, rho, phi, U, scheme);
```

其中：

- `rho`：面心密度场；
- `phi`：被输运变量；
- `U`：面心速度场；
- `scheme`：对流离散格式。

当前枚举格式：

```cpp
enum class DivType
{
    FUD,
    SUD,
    CD,
    MINMOD,
    MUSCL
};
```

当前主要实现：

| 格式 | 状态 | 说明 |
|---|---|---|
| `FUD` | 已实现 | 一阶迎风格式 |
| `SUD` | 已实现 | 二阶迎风格式，使用延迟修正 |
| `CD` | 已实现基础形式 | 中心差分格式 |
| `MINMOD` | 预留 | 尚未完整实现 |
| `MUSCL` | 预留 | 尚未完整实现 |

典型用法：

```cpp
FaceField<Scalar> rho("rho", &mesh);
FaceField<Vector<Scalar>> U("U", &mesh);
Field<Scalar> phi("T", &mesh);
SparseMatrix<Scalar> A(&mesh);

rho.setValue(1.0);
U.setValue(Vector<Scalar>(10, 0, 0));

fvm::Div(A, rho, phi, U, fvm::DivType::FUD);
```

---

### 2.3 源项 `fvm::Source`

源项接口分为两类。

#### 只含显式源项

```cpp
fvm::Source(matrix, S_c);
```

对应：

```math
S_\phi = S_c
```

#### 同时含线性源项和显式源项

```cpp
fvm::Source(matrix, S_p, S_c);
```

对应：

```math
S_\phi = S_p \phi + S_c
```

装配时：

- `S_c` 加入右端项；
- `S_p` 加入矩阵对角项。

典型用法：

```cpp
CellField<Scalar> Q("Q", &mesh);
Q.setValue(1000.0);

fvm::Source(A, Q);
```

---

## 3. 稀疏矩阵模块

稀疏矩阵类位于：

```text
src/core/Math/SparseMatrix.hpp
```

核心类型：

```cpp
template <typename Tp>
class SparseMatrix;
```

其中矩阵系数为 `Scalar`，右端项和未知量类型为模板参数 `Tp`。因此同一个矩阵结构可以服务于：

```cpp
SparseMatrix<Scalar>
SparseMatrix<Vector<Scalar>>
```

当前支持：

- 从普通二维矩阵构造；
- 从 `Mesh` 构造稀疏结构；
- CSR 压缩存储；
- 设置矩阵元素；
- 累加矩阵元素；
- 设置和累加右端项；
- 获取 CSR 三数组；
- 矩阵清零；
- 与 FVM 离散算子联动。

CSR 主要数组：

```text
values_      # 非零系数
colIndexs_   # 非零系数对应的列号
rowPointer_  # 每一行起始位置
```

示例：

```cpp
SparseMatrix<Scalar> A(&mesh);

A.setValue(0, 0, 10.0);
A.addValue(0, 1, -2.0);
A.addB(0, 5.0);

A.compress();
```

---

## 4. 线性求解器

线性求解器位于：

```text
src/core/Math/Solver.hpp
```

核心类型：

```cpp
template <typename Tp>
class Solver;
```

当前主要支持：

- Jacobi 迭代；
- 标量未知量；
- 矢量未知量；
- 串行求解；
- 基于线程池的并行 Jacobi 求解；
- 欠松弛；
- 最大迭代次数控制；
- 残差容差控制。

求解器枚举中预留了：

```cpp
GaussSeidel
AMG
```

但当前主要实现和使用的是：

```cpp
Solver<Tp>::Method::Jacobi
```

典型用法：

```cpp
Solver<Scalar> solver(A, Solver<Scalar>::Method::Jacobi, 100000);

solver.init(T.getCellField_0().getData());
solver.setTolerance(1e-10);
solver.relax(0.7);
solver.solve();
```

开启并行：

```cpp
solver.setParallel();
solver.solve();
```

并行求解当前通过 `ThreadPool` 完成，主要策略是将 CSR 矩阵行划分为多个连续区间，每个线程只写自己负责的 `x_[rowBegin, rowEnd)` 区间，从而减少多线程写冲突。

---

## 5. SIMPLE 算法

SIMPLE 算法位于：

```text
src/core/algorithm/SIMPLE/SIMPLE.hpp
```

命名空间：

```cpp
algorithm::simple
```

核心类型：

```cpp
algorithm::simple::SIMPLE
```

当前用于求解不可压缩稳态流动问题。基本未知量为：

- 速度场 `U`；
- 压力场 `p`。

物性参数为：

- 密度 `rho`；
- 粘度 `mu` 或运动粘度形式下的等效扩散系数。

---

### 5.1 SIMPLE 参数

主要参数位于：

```cpp
algorithm::simple::SIMPLE::Options
```

常用参数：

```cpp
algorithm::simple::SIMPLE::Options options;

options.maxOuterIterations = 5000;
options.convergenceTolerance = 1e-8;

options.momentumMaxIterations = 5000;
options.momentumTolerance = 1e-10;

options.pressureMaxIterations = 10000;
options.pressureTolerance = 1e-10;

options.alphaU = 0.7;
options.alphaP = 0.3;

options.nNonOrthogonalCorrectors = 2;
options.divScheme = fvm::DivType::SUD;

options.useRhieChow = true;
options.useParallel = true;
```

对于压力出口，可以设置：

```cpp
options.fixedPressurePatches = {"outlet"};
```

如果没有指定固定压力出口，程序会使用压力参考单元和惩罚项固定压力基准：

```cpp
options.pressureReferenceCell = 0;
options.pressureReferenceValue = 0.0;
options.pressureReferencePenalty = 1e20;
```

---

### 5.2 SIMPLE 计算流程

当前 SIMPLE 主流程大致为：

```text
1. U 和 p 从 cell 插值到 face
2. 计算 Rhie-Chow 面速度 Uf
3. 进入 SIMPLE 外迭代
4. 组装并求解动量方程
5. 计算 rAU = V / aP
6. 更新 Rhie-Chow 面速度
7. 组装并求解压力修正方程
8. 修正面速度
9. 修正压力
10. 修正速度
11. 计算连续性残差、速度残差、压力残差
12. 判断是否收敛
```

---

### 5.3 三维方腔 SIMPLE 示例

```cpp
#include <iostream>

#include "Field.hpp"
#include "Geometry/Mesh.h"
#include "SIMPLE.hpp"
#include "Vector.hpp"

using Scalar = double;

int main()
{
    Mesh mesh("tempFile/OpenFOAM_tutorials/cavity3D_4/constant/polyMesh");

    Field<Vector<Scalar>> U("U", &mesh);
    Field<Scalar> p("p", &mesh);

    U.setValue(Vector<Scalar>(0, 0, 0));
    p.setValue(0.0);

    U.setBoundaryCondition("movingWall", 1, 0, Vector<Scalar>(5, 0, 0));
    U.setBoundaryCondition("fixedWalls", 1, 0, Vector<Scalar>(0, 0, 0));

    p.setBoundaryCondition("movingWall", 0, 1, 0.0);
    p.setBoundaryCondition("fixedWalls", 0, 1, 0.0);

    FaceField<Scalar> rho("rho", &mesh);
    rho.setValue(1.0);

    FaceField<Scalar> nu("nu", &mesh);
    nu.setValue(0.01);

    algorithm::simple::SIMPLE::Options options;
    options.maxOuterIterations = 5000;
    options.alphaU = 0.7;
    options.alphaP = 0.3;
    options.convergenceTolerance = 1e-8;
    options.nNonOrthogonalCorrectors = 2;
    options.divScheme = fvm::DivType::SUD;
    options.useParallel = true;

    algorithm::simple::SIMPLE solver(U, p, rho, nu, options);
    solver.solve();

    U.writeToFile("U_SIMPLE.dat");
    p.writeToFile("p_SIMPLE.dat");

    return 0;
}
```

---

## 6. 项目目录结构

当前推荐目录结构如下：

```text
CFD_FVM_Solver
├── CMakeLists.txt
├── README.md
├── src
│   ├── CMakeLists.txt
│   ├── main.cpp
│   └── core
│       ├── CMakeLists.txt
│       ├── Geometry
│       │   ├── BoundaryPatch.cpp
│       │   ├── BoundaryPatch.h
│       │   ├── Cell.cpp
│       │   ├── Cell.h
│       │   ├── Face.cpp
│       │   ├── Face.h
│       │   ├── Mesh.cpp
│       │   └── Mesh.h
│       ├── Fields
│       │   ├── BaseField.hpp
│       │   ├── BoundaryCondition.hpp
│       │   ├── CellField.hpp
│       │   ├── FaceField.hpp
│       │   ├── Field.hpp
│       │   ├── FieldOperators
│       │   │   ├── Divergence.hpp
│       │   │   ├── Gradient.hpp
│       │   │   └── Laplacian.hpp
│       │   └── FieldType.hpp
│       ├── FVM
│       │   ├── Div.hpp
│       │   ├── DivType.h
│       │   ├── Laplacian.hpp
│       │   └── Source.hpp
│       ├── Math
│       │   ├── Solver.hpp
│       │   ├── SparseMatrix.hpp
│       │   ├── Tensor.hpp
│       │   ├── Vector.hpp
│       │   └── interpolation.hpp
│       ├── algorithm
│       │   └── SIMPLE
│       │       └── SIMPLE.hpp
│       └── thread_pool
│           └── threadpool.hpp
├── example
│   ├── CMakeLists.txt
│   ├── simple
│   │   ├── cavity2d_structured.cpp
│   │   ├── cavity2d_unstructured.cpp
│   │   ├── cavity3d_unstructured.cpp
│   │   └── pitz_daily_steady.cpp
│   ├── laplacian
│   │   ├── scalar_dirichlet.cpp
│   │   ├── scalar_neumann.cpp
│   │   └── vector_laplacian.cpp
│   ├── div
│   │   └── scalar_fud.cpp
│   ├── source
│   │   ├── heat_source.cpp
│   │   └── source_assembly_test.cpp
│   └── matrix
│       ├── sparse_matrix_access_test.cpp
│       └── sparse_matrix_solver_test.cpp
├── tempFile
├── build
├── bin
└── lib
```

---

## 7. 编译环境

### 7.1 基本要求

建议环境：

- Linux；
- GCC 11 或更新版本；
- CMake 3.15 或更新版本；
- C++20。

当前顶层 CMake 使用：

```cmake
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
```

---

### 7.2 编译

在项目根目录执行：

```bash
rm -rf build
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

Debug 编译：

```bash
rm -rf build
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build -j
```

---

### 7.3 输出目录

当前输出路径为：

```text
bin/    # 可执行文件
lib/    # 动态库 / 静态库
```

例如编译后可能生成：

```text
bin/CFD_FVM_Solver
bin/simple_cavity3d_unstructured
bin/simple_pitz_daily_steady
bin/simple_cavity2d_unstructured
bin/simple_cavity2d_structured
bin/laplacian_scalar_dirichlet
bin/laplacian_scalar_neumann
bin/laplacian_vector_laplacian
bin/div_scalar_fud
bin/source_heat_source
bin/source_source_assembly_test
bin/matrix_sparse_matrix_access_test
bin/matrix_sparse_matrix_solver_test
lib/libcore.so
```

---

## 8. 运行算例

### 8.1 main.cpp(测试专用)

```bash
./bin/CFD_FVM_Solver
```

### 8.2 运行三维非结构方腔算例

```bash
./bin/simple_cavity3d_unstructured
```

### 8.3 运行二维结构方腔算例

```bash
./bin/simple_cavity2d_structured
```

### 8.4 运行二维非结构方腔算例

```bash
./bin/simple_cavity2d_unstructured
```

### 8.5 运行台阶流算例

```bash
./bin/simple_pitz_daily_steady
```

### 8.6 运行扩散项测试

```bash
./bin/laplacian_scalar_dirichlet
./bin/laplacian_scalar_neumann
./bin/laplacian_vector_laplacian
```

### 8.7 运行对流项测试

```bash
./bin/div_scalar_fud
```

### 8.8 运行源项测试

```bash
./bin/source_heat_source
./bin/source_source_assembly_test
```

### 8.9 运行稀疏矩阵测试

```bash
./bin/matrix_sparse_matrix_access_test
./bin/matrix_sparse_matrix_solver_test
```

---

## 9. 添加新算例

项目中的 `example/CMakeLists.txt` 会自动递归收集 `example/` 下的所有 `.cpp` 文件。

因此新增算例只需要创建新的 `.cpp` 文件，例如：

```text
example/simple/my_new_case.cpp
```

然后重新配置和编译：

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

目标名会由相对路径自动生成：

```text
example/simple/my_new_case.cpp
```

会生成：

```text
bin/simple_my_new_case
```

---

## 10. 基本开发流程

### 10.1 创建网格

当前推荐使用 OpenFOAM 工具生成网格，例如：

```bash
blockMesh
```

或使用 Gmsh 生成 `.msh` 后转换：

```bash
gmshToFoam case.msh
```

然后在代码中读取：

```cpp
Mesh mesh("path/to/case/constant/polyMesh");
```

---

### 10.2 创建场

```cpp
Field<Scalar> p("p", &mesh);
Field<Vector<Scalar>> U("U", &mesh);

p.setValue(0.0);
U.setValue(Vector<Scalar>(0, 0, 0));
```

---

### 10.3 设置边界条件

```cpp
U.setBoundaryCondition("inlet", 1, 0, Vector<Scalar>(10, 0, 0));
U.setBoundaryCondition("outlet", 0, 1, Vector<Scalar>(0, 0, 0));
U.setBoundaryCondition("walls", 1, 0, Vector<Scalar>(0, 0, 0));

p.setBoundaryCondition("inlet", 0, 1, 0.0);
p.setBoundaryCondition("outlet", 1, 0, 0.0);
p.setBoundaryCondition("walls", 0, 1, 0.0);
```

---

### 10.4 设置物性参数

```cpp
FaceField<Scalar> rho("rho", &mesh);
rho.setValue(1.0);

FaceField<Scalar> nu("nu", &mesh);
nu.setValue(0.01);
```

---

### 10.5 调用 SIMPLE 求解器

```cpp
algorithm::simple::SIMPLE::Options options;

options.maxOuterIterations = 5000;
options.alphaU = 0.7;
options.alphaP = 0.3;
options.convergenceTolerance = 1e-8;
options.nNonOrthogonalCorrectors = 2;
options.divScheme = fvm::DivType::SUD;
options.useParallel = true;

algorithm::simple::SIMPLE solver(U, p, rho, nu, options);
solver.solve();
```

---

### 10.6 输出结果

```cpp
U.writeToFile("U_SIMPLE.dat");
p.writeToFile("p_SIMPLE.dat");
```

输出文件可使用 Tecplot 或兼容 Tecplot 数据格式的软件进行可视化。

---

## 11. 当前已包含的典型算例

### 11.1 SIMPLE 算例

| 文件 | 内容 |
|---|---|
| `example/simple/cavity3d_unstructured.cpp` | 三维非结构方腔 |
| `example/simple/cavity2d_unstructured.cpp` | 二维非结构方腔 |
| `example/simple/cavity2d_structured.cpp` | 二维结构方腔 |
| `example/simple/pitz_daily_steady.cpp` | 台阶流 / 后台阶流 |

### 11.2 Laplacian 算例

| 文件 | 内容 |
|---|---|
| `example/laplacian/scalar_dirichlet.cpp` | 标量扩散，第一类边界条件 |
| `example/laplacian/scalar_neumann.cpp` | 标量扩散，第二类边界条件 |
| `example/laplacian/vector_laplacian.cpp` | 矢量扩散测试 |

### 11.3 Div 算例

| 文件 | 内容 |
|---|---|
| `example/div/scalar_fud.cpp` | 标量对流项，一阶迎风格式 |

### 11.4 Source 算例

| 文件 | 内容 |
|---|---|
| `example/source/source_assembly_test.cpp` | 源项装配测试 |
| `example/source/heat_source.cpp` | 带源项导热问题 |

### 11.5 Matrix 算例

| 文件 | 内容 |
|---|---|
| `example/matrix/sparse_matrix_access_test.cpp` | 稀疏矩阵访问测试 |
| `example/matrix/sparse_matrix_solver_test.cpp` | 稀疏矩阵 Jacobi 求解测试 |

---

## 12. 代码设计说明

### 12.1 `core` 作为核心库

`src/core` 被编译为共享库：

```cmake
add_library(core SHARED ...)
```

上层主程序和所有 example 算例统一链接：

```cmake
target_link_libraries(target PRIVATE core)
```

这种结构的优点是：

- 算例代码和核心算法代码解耦；
- 新增算例不需要改动核心库；
- `src/main.cpp` 可以保留为主入口或开发调试入口；
- `example/` 可以专门管理不同验证算例。

---

### 12.2 模板化场类型

大量核心类采用模板：

```cpp
Field<Tp>
CellField<Tp>
FaceField<Tp>
SparseMatrix<Tp>
Solver<Tp>
```

因此同一套离散和求解流程可以处理：

```cpp
Scalar
Vector<Scalar>
Tensor<Scalar>
```

其中常用组合为：

```cpp
Field<Scalar>              // 压力、温度等标量场
Field<Vector<Scalar>>      // 速度等矢量场
SparseMatrix<Scalar>       // 标量方程
SparseMatrix<Vector<Scalar>> // 矢量方程
```

---

### 12.3 FVM 离散项与矩阵装配

FVM 算子不直接返回矩阵，而是向已有矩阵累加贡献：

```cpp
fvm::Div(A, rho, phi, U, scheme);
fvm::Laplacian(A, gamma, phi);
fvm::Source(A, S);
```

这使得多个项可以共同装配到同一个线性系统中：

```cpp
SparseMatrix<Scalar> A(&mesh);

fvm::Div(A, rho, T, U, fvm::DivType::FUD);
fvm::Laplacian(A, gamma, T);
fvm::Source(A, Q);
```

最终形成：

```math
A \phi = b
```

---

### 12.4 非正交修正

扩散项中将面向量分解为：

```math
\mathbf{S}_f = \mathbf{E}_f + \mathbf{T}_f
```

其中：

- `E_f`：正交贡献；
- `T_f`：非正交修正贡献。

当前实现中，非正交修正以显式源项形式加入右端项，用于处理非正交网格。

---

### 12.5 Rhie-Chow 插值

SIMPLE 算法中包含 Rhie-Chow 面速度修正，用于抑制 colocated 网格上速度-压力解耦导致的棋盘格压力振荡。

控制参数：

```cpp
options.useRhieChow = true;
```

当前 SIMPLE 流程中会计算面速度 `Uf_`，并在压力修正过程中修正面通量。

---

## 13. 当前限制

当前项目仍处于开发阶段，以下功能仍需完善：

1. **时间项离散尚未系统实现**

   当前主要面向稳态问题或伪时间迭代形式。

2. **线性求解器仍以 Jacobi 为主**

   `GaussSeidel` 和 `AMG` 在枚举中已有预留，但尚未完整实现。

3. **高阶对流格式仍需完善**

   `MINMOD` 和 `MUSCL` 已预留枚举，但当前尚未完整实现。

4. **湍流模型尚未实现**

   当前 SIMPLE 算例主要针对层流不可压缩流。

5. **边界条件类型仍需扩展**

   当前通过 Robin 统一形式覆盖常见 Dirichlet / Neumann / Robin 边界，复杂工程边界仍需继续封装。

6. **并行仅限共享内存线程并行**

   当前没有 MPI 区域分解并行。

---

## 14. 开发计划

后续计划包括：

- 完善 Gauss-Seidel / SOR 求解器；
- 引入 AMG 或其他更高效线性求解器；
- 增加瞬态项离散；
- 完善 `MINMOD / MUSCL` 等高分辨率对流格式；
- 完善压力出口、速度入口、壁面等边界条件封装；
- 增加湍流模型；
- 增加更多标准验证算例；
- 增加单元测试；
- 优化大规模网格下的内存和求解性能；
- 完善 Tecplot / VTK 等后处理输出格式。

---

## 15. Git 使用建议

新增或修改算例后，可以使用：

```bash
git status
git add .
git commit -m "Update examples and README"
```

推送到 GitHub：

```bash
git push
```

---

## 16. 联系方式

作者：XC
邮箱：yuanxch3@mail2.sysu.edu.cn

项目地址：

```text
https://github.com/yxcGi/CFD_FVM_Solver
```
