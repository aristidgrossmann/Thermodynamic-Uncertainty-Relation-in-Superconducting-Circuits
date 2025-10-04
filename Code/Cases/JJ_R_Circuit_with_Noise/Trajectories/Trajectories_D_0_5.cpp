#include <iostream>
#include <vector>           //dynamic array
#include <Eigen/Dense>      //for linear algebra matrix operations
#include <random>           //random generator
#include <functional>
#include <fstream>          // For file operations
#include <cmath>  
#include <vector>


#include "../../../SDE_Solver/include/SDE_Solver.hpp"

int main() {
    
    //////////////////////////////////  RSJ CIRCUIT FOKKER PLANCK MEAN & VARIANCE    /////////////////////////////////////////////////
    
    //Diffusion Coefficient
    const double D = 0.5;

    //Position Space discretization
    int N = 2e3;
    const double xmin = -20*2*pi;
    const double xmax = 30*2*pi;
    Vector x = Vector::LinSpaced(N, xmin, xmax);

    //initial PDF
    const double p_x_init = 0.0;
    int idx_init = (p_x_init-xmin)/(xmax-xmin)*N;
    Vector p0 = Vector::Zero(N);
    p0(idx_init) = 1.0/3.0/(x(1)-x(0));
    p0(idx_init-1) = 1.0/3.0/(x(1)-x(0));
    p0(idx_init+1) = 1.0/3.0/(x(1)-x(0));


    const int noof_timesteps = 1e4;
    const int update_frequency = 1e2;
    int noof_output_steps = noof_timesteps/update_frequency +1;

    double tolerance_CLT = 1e-6;
    double dt;

    //parameter study IV curve parameters
    const int noof_i_0_vals = 4;  //number of different parameter values
    double i_min = 1.25;
    double i_max = 2.0;
    Vector i_0_vals = Vector::LinSpaced(noof_i_0_vals, i_min, i_max);
    double i_0;

    //store Trajectories:
    Fokker_Planck_1D_result result;
    result.time = Vector(noof_output_steps); 
    result.mean = Vector(noof_output_steps);
    result.variance = Vector(noof_output_steps);

    for (int i = 0; i < noof_i_0_vals; ++i){

        printf ("completed: %d \n", i);
        i_0 = i_0_vals(i);
        std::cout << "i0 " << i_0<< std::endl;

        dt  = std::min((x(1)- x(0)) / (i_0 +1), std::pow((x(1)- x(0)), 2.0) / (2.0*D));

        Fokker_Planck_Drift_Function drift_function =  [i_0](double phi){
            return i_0 - std::sin(phi);
            // return i_0 ;
        };

        result = Fokker_Planck_FVM_1D_Crank_Nicolson_Solver(drift_function, D, x, p0, noof_timesteps, update_frequency, dt, tolerance_CLT);                       

        writeMatrixToFile(result.time, std::to_string(i).substr(0, 1) + "__time.txt");
        writeMatrixToFile(result.mean, std::to_string(i).substr(0, 1) + "__mean.txt");
        writeMatrixToFile(result.variance, std::to_string(i).substr(0, 1) + "__variance.txt");
        writeMatrixToFile(result.mean_velocity, std::to_string(i).substr(0, 1) + "__mean_velocity.txt");
        writeMatrixToFile(result.variance_velocity, std::to_string(i).substr(0, 1) + "__variance_velocity.txt");

    }
    return 0;
}
