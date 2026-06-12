# Changelog

## v0.1.0 - First runnable FVM solver release

### Added

- OpenFOAM-style `polyMesh` reading.
- Mesh topology and geometry calculation.
- Scalar, vector, and tensor field abstractions.
- Finite volume discretization operators:
  - `fvm::Laplacian`
  - `fvm::Div`
  - `fvm::Source`
- CSR sparse matrix assembly.
- Jacobi linear solver with scalar and vector unknown support.
- SIMPLE algorithm for steady incompressible flow.
- Rhie-Chow interpolation for collocated grids.
- Non-orthogonal correction for diffusion terms.
- Example cases:
  - 2D lid-driven cavity
  - 2D backward-facing step
  - 3D lid-driven cavity
  - Laplacian tests
  - convection tests
  - source term tests
  - sparse matrix tests

### Known limitations

- Mainly supports steady problems.
- Linear solver is still primarily Jacobi-based.
- Turbulence models are not implemented.
- Boundary conditions are still being extended.
- Parallelism is currently shared-memory based only.
