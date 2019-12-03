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

    //MatrixXd xs = MatrixXd::Random(10,2);
    MatrixXd xs(2,2);
    xs << 1,0,2,0;
    MatrixXd vs = MatrixXd::Zero(2,2);
    MatrixXd ms = VectorXd::Ones(2);

    SimualationParameters params;
    auto sim = VerletSimulation(params,xs,vs,ms);
    sim.step(1);
}
