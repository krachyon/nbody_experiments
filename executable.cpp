//
// Created by basti on 22.10.19.
//

#include <iostream>
#include "nbody.h"
#include <Eigen/LU>

int main()
{
    void threading();

    using namespace Eigen;

    MatrixXd xs = MatrixXd::Random(10,2);
    MatrixXd vs = MatrixXd::Zero(10,2);
    MatrixXd ms = VectorXd::Ones(10);

    std::cout << xs;
    std::cout << vs;
    std::cout << ms;

    step(xs,vs,ms,1);
}
