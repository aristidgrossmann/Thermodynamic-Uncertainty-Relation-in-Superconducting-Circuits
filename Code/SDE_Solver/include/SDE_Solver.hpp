// Libraries
#pragma once

#include <iostream>
#include <vector>           //dynamic array
#include <Eigen/Dense>      //for linear algebra matrix operations
#include <random>           //random generator
#include <functional>
#include <fstream>          // For file operations
#include <cmath>            // math
#include <vector>           // Vectors
#include <stdexcept>        // error handlind




constexpr double pi = 3.14159265358979323846;

// Define aliases
using namespace std;
using Vector = Eigen::VectorXd; 
using Matrix = Eigen::MatrixXd;
using DriftFunction = std::function<Matrix(const Matrix&)>;  // returns drift vectors (as rows of a matrix)
using DiffusionFunction = std::function<std::vector<Matrix>(const Matrix&)>; // returns diffusion matrices 
using Fokker_Planck_Drift_Function = std::function<double(double)>;

struct Cumulants{
    Matrix mean;
    Matrix mean_std;
    Matrix variance;
    Matrix variance_std;
};

struct Long_Time_Mean_Velocity{
    Matrix mean_velocity;
    Matrix mean_velocity_std;
};

struct SDE_1D_result{
    Vector time;
    Vector mean;
    Vector variance;
    Vector mean_velocity;
    Vector variance_velocity;
    Vector mean_std;
    Vector variance_std;
    Vector mean_velocity_std;
    Vector variance_velocity_std;
};

struct Fokker_Planck_1D_result{
    Vector time;
    Vector mean;
    Vector variance;
    Vector mean_velocity;
    Vector variance_velocity;
};

///////////////////////////   HELPER FUNCTIONS  ////////////////////////////////////
void writeMatrixToFile(const Matrix& matrix, const std::string& filename);

double integrate_simpson(const Vector& x, const Vector& y);

Vector thomas_algorithm(
    const Vector& sub_diagonal,  // Length N-1 
    const Vector& diagonal,      // Length N
    const Vector& super_diagonal,// Length N-1 
    const Vector& rhs            // Length N (RHS vector b)
);

//array<Matrix, noof_timesteps>
// Generate precomputed Brownian increments (shape: time_steps x M x num_samples)
std::vector<Matrix> generateWienerIncrements(const int& M, const int& noof_timesteps, const int& noof_MC_samples, const double& dt);




///////////////////////////   SDE SOLVERS (LANGEVIN)  ////////////////////////////////////

Cumulants SDE_EulerMaruyama_Solver(
    const int& N, const int& M, //N: state vector dimension, M: Wiener process dimension
    const Vector& X0,   //X0: initial state
    const DriftFunction& drift_batch,   
    const DiffusionFunction& diffusion_batch,
    const double& T, 
    const int& noof_timesteps, const int& noof_MC_samples);


void SDE_Trajectory_EulerMaruyama_Solver(
    const int& N, const int& M, //N: state vector dimension, M: Wiener process dimension
    const Vector& X0,   //X0: initial state
    const DriftFunction& drift_batch,   
    const DiffusionFunction& diffusion_batch,
    const double& T, 
    const int& noof_timesteps, const int& noof_MC_samples);


Cumulants SDE_Time_Averaged_EulerMaruyama_Solver(
    const int& N, const int& M, //N: state vector dimension, M: Wiener process dimension
    const Vector& X0,   //X0: initial state
    const DriftFunction& drift_batch,   
    const DiffusionFunction& diffusion_batch,
    const double& T, 
    const int& noof_timesteps, const int& noof_MC_samples);

Long_Time_Mean_Velocity SDE_Long_Time_EulerMaruyama_Solver(
    const int& N, const int& M, //N: state vector dimension, M: Wiener process dimension
    const Vector& X0,   //X0: initial state
    const DriftFunction& drift_batch,   
    const DiffusionFunction& diffusion_batch,
    const double& T, 
    const int& noof_timesteps);

SDE_1D_result SDE_1D_EulerMaruyama_Solver(
    Fokker_Planck_Drift_Function drift_function, 
    const double& D,
    const double& x0,   //X0: initial state 
    const int& noof_MC_samples,       //Number of samples 
    const int& noof_timesteps,
    const int& update_frequency,
    double dt, 
    double tolerance_CLT);




///////////////////////////   FOKKER PLANCK FRAMEWORK SOLVERS  ////////////////////////////////////
Fokker_Planck_1D_result Fokker_Planck_FVM_1D_Explicit_Euler_Solver(
    Fokker_Planck_Drift_Function drift_function,   
    const double D,                                 
    const Vector& x,                          
    const Vector& p0,                         
    const int& noof_timesteps,                
    const int& update_frequency,              
    double dt);

Fokker_Planck_1D_result Fokker_Planck_FVM_1D_Crank_Nicolson_Solver(
    Fokker_Planck_Drift_Function drift_function,
    double D,
    const Vector& x,
    const Vector& p0,
    const int& noof_timesteps,
    const int& update_frequency,
    double dt, 
    double tolerance_CLT);






