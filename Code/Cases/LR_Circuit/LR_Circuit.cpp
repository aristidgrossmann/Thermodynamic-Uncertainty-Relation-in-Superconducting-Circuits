#include <iostream>
#include <vector>           //dynamic array
#include <Eigen/Dense>      //for linear algebra matrix operations
#include <random>           //random generator
#include <functional>
#include <fstream>          // For file operations

#include "../../SDE_Solver/include/SDE_Solver.hpp"

int main() {
    
    //////////////////////////////////   LR CIRCUIT    /////////////////////////////////////////////////
    int N = 1;
    int M = 1;
    int noof_timesteps = 5e3;
    int noof_samples = 1e5;
    double T = 5.0;

    Vector X0 = Vector::Zero(N);

    // Vectorized drift: input N x samples, output N x samples
    DriftFunction drift = [](const Matrix& X) {
        int N = X.rows();
        int noof_samples = X.cols();
        
        Matrix result = -1.0*X;
        return result;
    };

    // Vectorized diffusion: returns a vector of M matrices (each N x samples)
    DiffusionFunction diffusion = [](const Matrix& X) {
        int N = X.rows();
        int noof_samples = X.cols();

        std::vector<Matrix> result(N);
        result[0] = Matrix::Ones(noof_samples, 1);   //generally: example:  result[n] = Matrix::Ones(noof_samples, M);
        return result;
    };

    SDE_Trajectory_EulerMaruyama_Solver(N, M, X0, drift, diffusion, T, noof_timesteps, noof_samples);
    

    
    return 0;
}
