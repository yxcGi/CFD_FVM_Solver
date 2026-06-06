#include <exception>
#include <iostream>
#include <vector>

#include "Solver.hpp"
#include "SparseMatrix.hpp"

using Scalar = double;

int main()
{
    try {
        std::vector<std::vector<Scalar>> A{
            {1, 1},
            {1, 2}
        };

        std::vector<Scalar> b{3, 5};

        SparseMatrix<Scalar> A_sparse(A);
        A_sparse.setB(b);
        A_sparse.printMatrix();
        A_sparse.compress();

        Solver<Scalar> solver(A_sparse, Solver<Scalar>::Method::Jacobi, 100);
        solver.init(b);
        solver.solve();

        A_sparse.printMatrix();
    }
    catch (const std::exception& e) {
        std::cerr << "Exception: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
