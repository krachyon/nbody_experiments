//
// Created by basti on 22.10.19.
//

#include <iostream>
#include "nbody.h"
#include <Eigen/LU>
#include <chrono>

using namespace std::chrono;
int main()
{
    void threading();

    using namespace Eigen;

    MatrixType xs = MatrixType::Random(200,2);
    //MatrixXd xs(2,2);
    //xs << 1,0,2,0;
    MatrixType vs = MatrixType::Zero(200,2);
    MatrixType ms = VectorType::Ones(500);


    MatrixType xs_b = MatrixType::Random(200,2);
    //MatrixXd xs(2,2);
    //xs << 1,0,2,0;
    MatrixType vs_b = MatrixType::Zero(200,2);
    MatrixType ms_b = VectorType::Ones(500);

    SimualationParameters params;

    auto leapfrog = LeapfrogSimulation(params,xs,vs,ms);
    auto euler = EulerSimulation(params,xs,vs,ms);

    auto start = high_resolution_clock::now();
    leapfrog.step(5000);
    auto end = high_resolution_clock::now();
    std::cout << "leapfrog: " << std::chrono::duration_cast<microseconds>(end-start).count() << std::endl;

    start = high_resolution_clock::now();
    euler.step(5000);
    end = high_resolution_clock::now();
    std::cout << "euler: " << std::chrono::duration_cast<microseconds>(end-start).count() << std::endl;

    start = high_resolution_clock::now();
    leapfrog.step(5000);
    end = high_resolution_clock::now();
    std::cout << "leapfrog: " << std::chrono::duration_cast<microseconds>(end-start).count() << std::endl;

    start = high_resolution_clock::now();
    euler.step(5000);
    end = high_resolution_clock::now();
    std::cout << "euler: " << std::chrono::duration_cast<microseconds>(end-start).count() << std::endl;

}
