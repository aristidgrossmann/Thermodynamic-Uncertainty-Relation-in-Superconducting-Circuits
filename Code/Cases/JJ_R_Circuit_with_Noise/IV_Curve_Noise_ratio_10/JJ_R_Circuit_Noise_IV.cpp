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

    //Simulation parameters
    int N = 1;
    int M = 1;
    double T = 1e7;
    int noof_timesteps = 1e6;

    //parameter study IV curve parameters
    int noof_relative_current_values = 201;  //number of different parameter values
    double i_min = 0.0;
    double i_max = 2.0;
    Vector relative_current_values = Vector::LinSpaced(noof_relative_current_values, i_min, i_max);
    constexpr double pi = 3.14159265358979323846;
    double relative_current;
    double noise_ratio = 10;

    //initial condition
    Vector X0 = Vector::Zero(N);


    //store I-V curve:
    Cumulants IV_curve;
    IV_curve.mean = Matrix::Ones(N, noof_relative_current_values);
    IV_curve.mean_std = Matrix::Ones(N, noof_relative_current_values);
    IV_curve.variance = Matrix::Ones(N, noof_relative_current_values);
    IV_curve.variance_std = Matrix::Ones(N, noof_relative_current_values);

    // Vectorized diffusion: returns a vector of M matrices (each N x samples)
    DiffusionFunction diffusion = [noise_ratio](const Matrix& X) {
        int N = X.rows();
        int noof_samples = X.cols();

        std::vector<Matrix> result(N);
        result[0] = Matrix::Constant(noof_samples, 1, noise_ratio);   //generally: example:  result[n] = Matrix::Ones(noof_samples, M);
        return result;
    };

    for (int i = 0; i < noof_relative_current_values; ++i){

        printf ("completed: %d \n", i);
        relative_current = relative_current_values(i);

        DriftFunction drift = [relative_current](const Matrix& X) {
            int N = X.rows();
            int noof_samples = X.cols();
            
            Matrix result = Matrix::Constant(N, noof_samples, relative_current) - X.array().sin().matrix();
            return result;
        };

        Long_Time_Mean_Velocity Long_Time_results = SDE_Long_Time_EulerMaruyama_Solver(
            N, M, X0, drift, diffusion, T, noof_timesteps);

        IV_curve.mean.col(i) = Long_Time_results.mean_velocity;
        IV_curve.mean_std.col(i) = Long_Time_results.mean_velocity_std;
    }
    
    writeMatrixToFile(relative_current_values, "Relative_Current_IV.txt");
    writeMatrixToFile(IV_curve.mean.transpose(), "Voltage_mean_IV.txt");
    writeMatrixToFile(IV_curve.mean_std.transpose(), "Voltage_mean_std_IV.txt");
    return 0;
}
