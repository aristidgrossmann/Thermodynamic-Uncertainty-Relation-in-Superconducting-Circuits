// Libraries
#pragma once

#include <iostream>
#include <vector>           //dynamic array
#include <Eigen/Dense>      //for linear algebra matrix operations
#include <random>           //random generator
#include <functional>
#include <fstream>          // For file operations


// Define aliases
using Vector = Eigen::VectorXd; 
using Matrix = Eigen::MatrixXd;
using DriftFunction = std::function<Matrix(const Matrix&)>;  // returns drift vectors (as rows of a matrix)
using DiffusionFunction = std::function<std::vector<Matrix>(const Matrix&)>; // returns diffusion matrices 

using Cumulants = struct {
    Matrix mean;
    Matrix mean_std;
    Matrix variance;
    Matrix variance_std;
};


//array<Matrix, noof_timesteps>
// Generate precomputed Brownian increments (shape: time_steps x M x num_samples)
std::vector<Matrix> generateWienerIncrements(int M, int noof_timesteps, int noof_samples, double dt);


Cumulants SDE_EulerMaruyama_Solver(
    int N, int M, //N: state vector dimension, M: Wiener process dimension
    const Vector& X0,   //X0: initial state
    DriftFunction drift_batch,   
    DiffusionFunction diffusion_batch,
    double T, 
    int noof_timesteps, int noof_samples);


void writeMatrixToFile(const Eigen::MatrixXd& matrix, const std::string& filename);

void SDE_Trajectory_EulerMaruyama_Solver(
    int N, int M, //N: state vector dimension, M: Wiener process dimension
    const Vector& X0,   //X0: initial state
    DriftFunction drift_batch,   
    DiffusionFunction diffusion_batch,
    double T, 
    int noof_timesteps, int noof_samples);


Cumulants SDE_Time_Averaged_EulerMaruyama_Solver(
    int N, int M, //N: state vector dimension, M: Wiener process dimension
    const Vector& X0,   //X0: initial state
    DriftFunction drift_batch,   
    DiffusionFunction diffusion_batch,
    double T, 
    int noof_timesteps, int noof_samples);

