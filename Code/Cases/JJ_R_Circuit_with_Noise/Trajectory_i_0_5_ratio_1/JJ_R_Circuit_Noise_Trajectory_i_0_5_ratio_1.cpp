#include <iostream>
#include <vector>           //dynamic array
#include <Eigen/Dense>      //for linear algebra matrix operations
#include <random>           //random generator
#include <functional>
#include <fstream>          // For file operations
#include <cmath>  


#include "../../../SDE_Solver/include/SDE_Solver.hpp"

int main() {
    
    //////////////////////////////////  JJ R CIRCUIT    /////////////////////////////////////////////////
    constexpr double pi = 3.14159265358979323846;

    double relative_current = 0.5;
    double noise_ratio = 0.27;

    //Simulation parameters
    int N = 1;
    int M = 1;
    int noof_timesteps = 1e4;
    int noof_samples = 1e4;   //deterministic, so only 1 sample
    double T = 1e4;


    //initial condition
    Vector X0 = Vector::Constant(N, 0.0);


    // Vectorized drift: input N x samples, output N x samples
    DriftFunction drift = [relative_current](const Matrix& X) {
        int N = X.rows();
        int noof_samples = X.cols();
        
        Matrix result = Matrix::Constant(N, noof_samples, relative_current) - X.array().sin().matrix();
        return result;
    };

    // Vectorized diffusion: returns a vector of M matrices (each N x samples)
    DiffusionFunction diffusion = [noise_ratio](const Matrix& X) {
        int N = X.rows();
        int noof_samples = X.cols();

        std::vector<Matrix> result(N);
        result[0] = Matrix::Constant(noof_samples, 1, noise_ratio);   //generally: example:  result[n] = Matrix::Ones(noof_samples, M);
        return result;
    };

    SDE_Trajectory_EulerMaruyama_Solver(N, M, X0, drift, diffusion, T, noof_timesteps, noof_samples);

    return 0;
}
