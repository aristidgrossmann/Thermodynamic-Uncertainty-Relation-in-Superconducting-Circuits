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
    
    // Currents 
    const int noof_i_0_vals = 20;  
    double i_min = 0.1;
    double i_max = 2.0;
    Vector i_0_vals = Vector::LinSpaced(noof_i_0_vals, i_min, i_max);
    writeMatrixToFile(i_0_vals,  "i_0_vals_numerical.txt");
    double i_0;

    //Diffusion Coefficient
    Vector D_vals{{5, 10}}; 
    double D;

    //Position Space discretization
    int N = 8e3;
    const double xmin = -40*2*pi;
    const double xmax = 160*2*pi;
    Vector x = Vector::LinSpaced(N, xmin, xmax);

    //initial PDF
    const double p_x_init = 0.0;
    int idx_init = (p_x_init-xmin)/(xmax-xmin)*N;
    Vector p0 = Vector::Zero(N);
    p0(idx_init) = 1.0/3.0/(x(1)-x(0));
    // p0(idx_init-1) = 1.0/3.0/(x(1)-x(0));
    // p0(idx_init+1) = 1.0/3.0/(x(1)-x(0));


    const int noof_timesteps = 1e4;
    const int update_frequency = 1e2;
    int noof_output_steps = noof_timesteps/update_frequency +1;

    double dt;
    double tolerance_CLT = 1e-6;


    // store mean_velocity and variance_velocity for all D and i_0:
    Vector mean_velocities = Vector::Zero(noof_i_0_vals);
    Vector variance_velocities = Vector::Zero(noof_i_0_vals);
    //store Trajectories:
    Fokker_Planck_1D_result result;
    result.time = Vector(noof_output_steps); 
    result.mean = Vector(noof_output_steps);
    result.variance = Vector(noof_output_steps);

    for (int j = 0; j < D_vals.size(); ++j){
        D = D_vals(j);

        for (int i = 0; i < noof_i_0_vals; ++i){
            i_0 = i_0_vals(i);
            
            dt  = std::min((x(1)- x(0)) / (i_0 +1.0), std::pow((x(1)- x(0)), 2.0) / (2.0*D));

            Fokker_Planck_Drift_Function drift_function =  [i_0](double phi){
                return i_0 - std::sin(phi);
                // return i_0 ;
            };

            result = Fokker_Planck_FVM_1D_Crank_Nicolson_Solver(drift_function, D, x, p0, noof_timesteps, update_frequency, dt, tolerance_CLT);  
            
            mean_velocities(i) = result.mean_velocity(result.mean_velocity.size()-1);
            variance_velocities(i) = result.variance_velocity(result.variance_velocity.size()-1);

            printf ("completed: %d \n", i+j*noof_i_0_vals);
        }

        writeMatrixToFile(mean_velocities, std::to_string(D).substr(0, 3) + "_mean_velocities.txt");
        writeMatrixToFile(variance_velocities, std::to_string(D).substr(0, 3) + "_variance_velocities.txt");
    }
    return 0;
}
