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
    int noof_timesteps = 1e4;
    int noof_samples = 1;   //deterministic, so only 1 sample
    double T;

    //parameter study IV curve parameters
    int noof_relative_current_values = 1000;  //number of different parameter values
    int i_min = 0.0;
    int i_max = 5.0;
    Vector relative_current_values = Vector::LinSpaced(noof_relative_current_values, i_min, i_max);
    constexpr double pi = 3.14159265358979323846;
    double relative_current;

    //initial condition
    Vector X0 = Vector::Zero(N);


    //store I-V curve:
    Cumulants IV_curve;
    IV_curve.mean = Matrix::Ones(N, noof_relative_current_values);
    IV_curve.mean_std = Matrix::Ones(N, noof_relative_current_values);
    IV_curve.variance = Matrix::Ones(N, noof_relative_current_values);
    IV_curve.variance_std = Matrix::Ones(N, noof_relative_current_values);

    // Vectorized diffusion: returns a vector of M matrices (each N x samples)
    DiffusionFunction diffusion = [](const Matrix& X) {
        int N = X.rows();
        int noof_samples = X.cols();

        std::vector<Matrix> result(N);
        result[0] = Matrix::Zero(noof_samples, 1);   //generally: example:  result[n] = Matrix::Ones(noof_samples, M);
        return result;
    };

    for (int i = 0; i < noof_relative_current_values; ++i){

        relative_current = relative_current_values(i);

        if (relative_current <= 1.0){
            T = 1e3;
            X0 = Vector::Constant(N, std::asin(relative_current));
        }
        else{
            T = 2*pi/std::sqrt(std::pow(relative_current, 2.0) -1.0);
            X0 = Vector::Constant(N, -pi);
            
        }

        // Vectorized drift: input N x samples, output N x samples
        DriftFunction drift = [relative_current](const Matrix& X) {
            int N = X.rows();
            int noof_samples = X.cols();
            
            Matrix result = Matrix::Constant(N, noof_samples, relative_current) - X.array().sin().matrix();
            return result;
        };

        Cumulants Trajectory_results = SDE_EulerMaruyama_Solver(
            N, M, X0, drift, diffusion, T, noof_timesteps, noof_samples);

        IV_curve.mean.col(i) = (Trajectory_results.mean.col(Trajectory_results.mean.cols()-1).array()
                                - Trajectory_results.mean.col(0).array())/T;
        IV_curve.mean_std.col(i) = (Trajectory_results.mean_std.col(Trajectory_results.mean_std.cols()-1).array().square() 
                                    + Trajectory_results.mean_std.col(0).array().square()).sqrt()/T;
    }
    
    writeMatrixToFile(relative_current_values, "Relative_Current_IV.txt");
    writeMatrixToFile(IV_curve.mean.transpose(), "Voltage_mean_IV.txt");
    writeMatrixToFile(IV_curve.mean_std.transpose(), "Voltage_mean_std_IV.txt");
    return 0;
}
