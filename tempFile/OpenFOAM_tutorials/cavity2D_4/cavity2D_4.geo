// ================================================================
// 2D square cavity mesh for OpenFOAM
//
// Domain:
//   x in [0, 1]
//   y in [0, 1]
//   z thickness = 0.01
//
// Mesh:
//   20 x 20 structured square cells in x-y plane
//   extruded to one layer in z direction for OpenFOAM
//
// Boundary names:
//   top
//   bottom
//   left
//   rightWalls
//   frontAndBack
//
// Usage:
//   gmsh -3 cavity2D_quad_20x20.geo -format msh2 -o cavity2D_quad_20x20.msh
//   gmshToFoam cavity2D_quad_20x20.msh
// ================================================================

SetFactory("Built-in");

// ------------------------------------------------
// Basic parameters
// ------------------------------------------------
L = 1.0;
thickness = 0.01;

// 20 x 20 cells means 21 points on each edge
n = 21;

// OpenFOAM gmshToFoam is more stable with msh2 format
Mesh.MshFileVersion = 2.2;

// ------------------------------------------------
// Points
// ------------------------------------------------
Point(1) = {0.0, 0.0, 0.0};
Point(2) = {L,   0.0, 0.0};
Point(3) = {L,   L,   0.0};
Point(4) = {0.0, L,   0.0};

// ------------------------------------------------
// Lines
// ------------------------------------------------
Line(1) = {1, 2};   // bottom
Line(2) = {2, 3};   // right
Line(3) = {3, 4};   // top
Line(4) = {4, 1};   // left

// ------------------------------------------------
// Surface
// ------------------------------------------------
Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// ------------------------------------------------
// Structured 20 x 20 quadrilateral mesh
// ------------------------------------------------
Transfinite Curve {1, 3} = n;
Transfinite Curve {2, 4} = n;

Transfinite Surface {1} = {1, 2, 3, 4};

// Recombine triangles into quadrilaterals
Recombine Surface {1};

// ------------------------------------------------
// Extrude to one-cell-thick 3D mesh for OpenFOAM
//
// ext[0] : back surface, z = thickness
// ext[1] : volume
// ext[2] : surface generated from Line(1), bottom
// ext[3] : surface generated from Line(2), right
// ext[4] : surface generated from Line(3), top
// ext[5] : surface generated from Line(4), left
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

// z-normal patches for 2D calculation
Physical Surface("frontAndBack") = {1, ext[0]};

// fluid volume
Physical Volume("fluid") = {ext[1]};

// Generate 3D mesh
Mesh 3;