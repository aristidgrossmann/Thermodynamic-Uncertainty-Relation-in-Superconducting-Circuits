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
    
    //////////////////////////////////  RSJ CIRCUIT SDE MEAN & VARIANCE    /////////////////////////////////////////////////
    
    // Currents 
    const int noof_i_0_vals = 20;  
    double i_min = 0.1;
    double i_max = 2.0;
    Vector i_0_vals = Vector::LinSpaced(noof_i_0_vals, i_min, i_max);
    writeMatrixToFile(i_0_vals,  "i_0_vals_numerical.txt");
    double i_0;

    //Diffusion Coefficient
    Vector D_vals{{10.0}}; 
    double D;

    int noof_MC_samples = 5e4;
    const double x0 = 0.0; //initial state


    const int noof_timesteps = 1e4;
    const int update_frequency = 1e3;

    double dt;
    double tolerance_CLT = 1e-2;


    // store mean_velocity and variance_velocity for all D and i_0:
    Vector mean_velocities = Vector::Zero(noof_i_0_vals);
    Vector mean_velocities_std = Vector::Zero(noof_i_0_vals);

    Vector variance_velocities = Vector::Zero(noof_i_0_vals);
    Vector variance_velocities_std = Vector::Zero(noof_i_0_vals);

    //store Trajectories:
    SDE_1D_result result;

    for (int j = 0; j < D_vals.size(); ++j){
        D = D_vals(j);

        for (int i = 0; i < noof_i_0_vals; ++i){
            i_0 = i_0_vals(i);
            
            // dt  = std::min(1/ (i_0 +1.0), std::pow((0.1), 2.0) / (2.0*D));
            dt = 0.1/(i_0+1.0);
            std::cout << "Time step size: " << dt << std::endl;

            Fokker_Planck_Drift_Function drift_function =  [i_0](double phi){
                return i_0 - std::sin(phi);
                // return i_0 ;
            };

            result = SDE_1D_EulerMaruyama_Solver(drift_function, D, x0, noof_MC_samples, noof_timesteps, update_frequency, dt, tolerance_CLT);  
            
            mean_velocities(i) = result.mean_velocity(result.mean_velocity.size()-1);
            variance_velocities(i) = result.variance_velocity(result.variance_velocity.size()-1);

            mean_velocities_std(i) = result.mean_velocity_std(result.mean_velocity_std.size()-1);
            variance_velocities_std(i) = result.variance_velocity_std(result.variance_velocity_std.size()-1);

            printf ("completed: %d \n", i+j*noof_i_0_vals);
        }

        writeMatrixToFile(mean_velocities, std::to_string(D).substr(0, 3) + "_mean_velocities.txt");
        writeMatrixToFile(variance_velocities, std::to_string(D).substr(0, 3) + "_variance_velocities.txt");
        writeMatrixToFile(mean_velocities_std, std::to_string(D).substr(0, 3) + "_mean_velocities_std.txt");
        writeMatrixToFile(variance_velocities_std, std::to_string(D).substr(0, 3) + "_variance_velocities_std.txt");
    }
    return 0;
}
