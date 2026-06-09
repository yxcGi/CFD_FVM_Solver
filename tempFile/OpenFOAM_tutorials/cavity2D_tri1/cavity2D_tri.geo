// ================================================================
// 2D square cavity triangular mesh for OpenFOAM
//
// Domain:
//   x in [0, 1]
//   y in [0, 1]
//   thickness in z direction = 0.01
//
// Boundary names:
//   top
//   bottom
//   left
//   rightWalls
//   frontAndBack
//
// Usage:
//   gmsh -3 cavity2D_tri.geo -format msh2 -o cavity2D_tri.msh
//   gmshToFoam cavity2D_tri.msh
// ================================================================

SetFactory("Built-in");

// ------------------------------------------------
// Mesh parameters
// ------------------------------------------------
lc = 0.05;          // 网格尺寸，越小网格越密
thickness = 0.01;   // z方向厚度，一层网格，用于二维计算

// OpenFOAM 的 gmshToFoam 对 msh2 格式兼容性更好
Mesh.MshFileVersion = 2.2;

// 使用 Delaunay 三角网格算法
Mesh.Algorithm = 5;

// 保持三角形网格，不重组为四边形
Mesh.RecombineAll = 0;

// ------------------------------------------------
// Points: square [0,1] x [0,1]
// ------------------------------------------------
Point(1) = {0.0, 0.0, 0.0, lc};
Point(2) = {1.0, 0.0, 0.0, lc};
Point(3) = {1.0, 1.0, 0.0, lc};
Point(4) = {0.0, 1.0, 0.0, lc};

// ------------------------------------------------
// Boundary lines
// ------------------------------------------------
Line(1) = {1, 2};   // bottom: y = 0
Line(2) = {2, 3};   // right : x = 1
Line(3) = {3, 4};   // top   : y = 1
Line(4) = {4, 1};   // left  : x = 0

// ------------------------------------------------
// 2D surface
// ------------------------------------------------
Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// ------------------------------------------------
// Extrude 2D triangular mesh to one-layer 3D mesh
//
// ext[0]: top surface after extrusion, z = thickness
// ext[1]: volume
// ext[2]: side surface generated from Line(1), bottom
// ext[3]: side surface generated from Line(2), right
// ext[4]: side surface generated from Line(3), top
// ext[5]: side surface generated from Line(4), left
//
// Recombine here means:
//   triangle -> triangular prism
// ------------------------------------------------
ext[] = Extrude {0, 0, thickness}
{
    Surface{1};
    Layers{1};
    Recombine;
};

// ------------------------------------------------
// Physical groups for OpenFOAM patches
// ------------------------------------------------
Physical Surface("bottom") = {ext[2]};
Physical Surface("rightWalls") = {ext[3]};
Physical Surface("top") = {ext[4]};
Physical Surface("left") = {ext[5]};

// z = 0 and z = thickness surfaces
Physical Surface("frontAndBack") = {1, ext[0]};

// Fluid volume
Physical Volume("fluid") = {ext[1]};

// Generate 3D mesh
Mesh 3;