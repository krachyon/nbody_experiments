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

    MatrixXd xs = MatrixXd::Random(100,2);
    //MatrixXd xs(2,2);
    //xs << 1,0,2,0;
    MatrixXd vs = MatrixXd::Zero(100,2);
    MatrixXd ms = VectorXd::Ones(100);

    SimualationParameters params;
    //auto sim = VerletSimulation(params,xs,vs,ms);
    //sim.step(1);

    auto leapfrog = LeapfrogSimulation(params,xs,vs,ms);
    leapfrog.step(5000);
    std::cout << xs;
}
