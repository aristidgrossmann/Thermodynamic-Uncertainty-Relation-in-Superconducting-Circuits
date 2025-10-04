// Libraries
#include "../include/SDE_Solver.hpp"


//array<Matrix, noof_timesteps>
// Generate precomputed Brownian increments (shape: time_steps x M x num_samples)
std::vector<Matrix> generateWienerIncrements(int M, int noof_timesteps, int noof_samples, double dt) {
    std::mt19937 random_number_generator(std::random_device{}());  //define random number generator object of type mt19937
    std::normal_distribution<> normal_distribution(0.0, 1.0);      //define normal distribution object

    std::vector<Matrix> dWs(noof_timesteps); //generate noof_timesteps matrices with dimension (M x noof_samples)
    for (int t = 0; t < noof_timesteps; ++t) {    //loop over timesteps
        Matrix dW = Matrix::Zero(M, noof_samples);          //dW: contains Wiener increments (M x noof_samples)
        for (int m = 0; m < M; ++m) {             // loop over Wiener dimensions
            for (int s = 0; s < noof_samples; ++s) {   // loop over samples
                dW(m, s) = normal_distribution(random_number_generator);  //calculate Wiener increment
            }
        }
        dWs[t] = dW* std::sqrt(dt);
    }
    return dWs;
}


Cumulants SDE_EulerMaruyama_Solver(
    int N, int M, //N: state vector dimension, M: Wiener process dimension
    const Vector& X0,   //X0: initial state
    DriftFunction drift_batch,   
    DiffusionFunction diffusion_batch,
    double T, 
    int noof_timesteps, int noof_samples) {

    Cumulants result;
    result.mean = Matrix::Ones(N, noof_timesteps);  //store mean
    result.variance = Matrix::Ones(N, noof_timesteps); //store variance
    result.mean_std = Matrix::Ones(N, noof_timesteps);  //store std of the mean
    result.variance_std = Matrix::Ones(N, noof_timesteps); //store std of the variance

    Matrix X_centered = Matrix::Ones(N, noof_timesteps);  //store centered mean
    Matrix curtosis = Matrix::Ones(N, noof_timesteps); //store curtosis

    double dt = T / noof_timesteps;

    std::vector<Matrix> dWs = generateWienerIncrements(M, noof_timesteps, noof_samples, dt);  //pre-compute Wiener increments

    Matrix X = X0.replicate(1, noof_samples);  //matrix of initial states

    for (int t = 0; t < noof_timesteps; ++t) {  //loop over noof_timesteps
        
        //Write means first before update(X0)
        result.mean.col(t) = X.rowwise().mean();
        X_centered = X.colwise() - result.mean.col(t);
        //variance.col(t) = (X.colwise() - mean.col(t)).array().square().rowwise().sum()/(noof_samples-1);
        result.variance.col(t) = X_centered.array().square().rowwise().sum() / (noof_samples - 1);

        result.mean_std.col(t) = (result.variance.col(t)/noof_samples).array().sqrt();
        curtosis.col(t) = X_centered.array().pow(4).rowwise().sum() / noof_samples;
        result.variance_std.col(t) = ((curtosis.col(t).array() - result.variance.col(t).array().square()) / noof_samples).sqrt();
        //variance_std.col(t) = ((curtosis.col(t).colwise() - variance.col(t).array().square().matrix())/noof_samples).array().sqrt();

    
        double time = t * dt; //update time

        Matrix drift = drift_batch(X) * dt;  //compute drift for all samples (noof_samples)

        std::vector<Matrix> sigma = diffusion_batch(X); //compute diffusion matrix(size M, each N x noof_samples)

        // Compute diffusion contribution: N x num_samples
        Matrix diffusion = Matrix::Zero(N, noof_samples);
        for (int n = 0; n < N; ++n) {
            for (int s = 0; s < noof_samples; ++s){
                diffusion(n,s) = (sigma[n].row(s)).dot(dWs[t].col(s));
            }
        }

        X += drift + diffusion;   //update X (matrix N x noof_samples)
    }
    return result;
}


void writeMatrixToFile(const Eigen::MatrixXd& matrix, const std::string& filename) {
    // Open file in truncation mode (overwrites existing file or creates new)

    std::ofstream file(filename, std::ios::trunc);
    
    if (file.is_open()) {
        file << matrix;  // Eigen's built-in output formatting
        std::cout << "Successfully wrote to " << filename << "\n";
    } else {
        std::cerr << "Error: Failed to create/write to " << filename << "\n";
    }
}

// ===== Vectorized Monte Carlo Simulation =====
void SDE_Trajectory_EulerMaruyama_Solver(
    int N, int M, //N: state vector dimension, M: Wiener process dimension
    const Vector& X0,   //X0: initial state
    DriftFunction drift_batch,   
    DiffusionFunction diffusion_batch,
    double T, 
    int noof_timesteps, int noof_samples) {

    Matrix mean = Matrix::Ones(N, noof_timesteps);  //store mean
    Matrix X_centered = Matrix::Ones(N, noof_timesteps);  //store centered mean
    Matrix variance = Matrix::Ones(N, noof_timesteps); //store variance
    Matrix mean_std = Matrix::Ones(N, noof_timesteps);  //store std of the mean
    Matrix variance_std = Matrix::Ones(N, noof_timesteps); //store std of the variance
    Matrix curtosis = Matrix::Ones(N, noof_timesteps); //store std of the variance

    double dt = T / noof_timesteps;

    std::vector<Matrix> dWs = generateWienerIncrements(M, noof_timesteps, noof_samples, dt);  //pre-compute Wiener increments

    Matrix X = X0.replicate(1, noof_samples);  //matrix of initial states

    for (int t = 0; t < noof_timesteps; ++t) {  //loop over noof_timesteps
        
        //Write means first before update(X0)
        mean.col(t) = X.rowwise().mean();
        X_centered = X.colwise() - mean.col(t);
        //variance.col(t) = (X.colwise() - mean.col(t)).array().square().rowwise().sum()/(noof_samples-1);
        variance.col(t) = X_centered.array().square().rowwise().sum() / (noof_samples - 1);

        mean_std.col(t) = (variance.col(t)/noof_samples).array().sqrt();
        curtosis.col(t) = X_centered.array().pow(4).rowwise().sum() / noof_samples;
        variance_std.col(t) = ((curtosis.col(t).array() - variance.col(t).array().square()) / noof_samples).sqrt();
        //variance_std.col(t) = ((curtosis.col(t).colwise() - variance.col(t).array().square().matrix())/noof_samples).array().sqrt();

    
        double time = t * dt; //update time

        Matrix drift = drift_batch(X) * dt;  //compute drift for all samples (noof_samples)

        std::vector<Matrix> sigma = diffusion_batch(X); //compute diffusion matrix(size M, each N x noof_samples)

        // Compute diffusion contribution: N x num_samples
        Matrix diffusion = Matrix::Zero(N, noof_samples);
        for (int n = 0; n < N; ++n) {
            for (int s = 0; s < noof_samples; ++s){
                diffusion(n,s) = (sigma[n].row(s)).dot(dWs[t].col(s));
            }
        }

        X += drift + diffusion;   //update X (matrix N x noof_samples)
    }


    writeMatrixToFile(Vector::LinSpaced(noof_timesteps, 0, T), "time.txt");
    writeMatrixToFile(mean.transpose(), "mean.txt");
    writeMatrixToFile(mean_std.transpose(), "mean_std.txt");
    writeMatrixToFile(variance.transpose(), "variance.txt");
    writeMatrixToFile(variance_std.transpose(), "variance_std.txt");
}


Cumulants SDE_Time_Averaged_EulerMaruyama_Solver(
    int N, int M, //N: state vector dimension, M: Wiener process dimension
    const Vector& X0,   //X0: initial state
    DriftFunction drift_batch,   
    DiffusionFunction diffusion_batch,
    double T, 
    int noof_timesteps, int noof_samples) {

    Matrix mean = Matrix::Ones(N, noof_timesteps);  //store mean
    Matrix X_centered = Matrix::Ones(N, noof_timesteps);  //store centered mean
    Matrix variance = Matrix::Ones(N, noof_timesteps); //store variance
    Matrix mean_std = Matrix::Ones(N, noof_timesteps);  //store std of the mean
    Matrix variance_std = Matrix::Ones(N, noof_timesteps); //store std of the variance
    Matrix curtosis = Matrix::Ones(N, noof_timesteps); //store std of the variance

    double dt = T / noof_timesteps;

    std::vector<Matrix> dWs = generateWienerIncrements(M, noof_timesteps, noof_samples, dt);  //pre-compute Wiener increments

    Matrix X = X0.replicate(1, noof_samples);  //matrix of initial states

    for (int t = 0; t < noof_timesteps; ++t) {  //loop over noof_timesteps
        
        //Write means first before update(X0)
        mean.col(t) = X.rowwise().mean();
        X_centered = X.colwise() - mean.col(t);
        //variance.col(t) = (X.colwise() - mean.col(t)).array().square().rowwise().sum()/(noof_samples-1);
        variance.col(t) = X_centered.array().square().rowwise().sum() / (noof_samples - 1);

        mean_std.col(t) = (variance.col(t)/noof_samples).array().sqrt();
        curtosis.col(t) = X_centered.array().pow(4).rowwise().sum() / noof_samples;
        variance_std.col(t) = ((curtosis.col(t).array() - variance.col(t).array().square()) / noof_samples).sqrt();
        //variance_std.col(t) = ((curtosis.col(t).colwise() - variance.col(t).array().square().matrix())/noof_samples).array().sqrt();

    
        double time = t * dt; //update time

        Matrix drift = drift_batch(X) * dt;  //compute drift for all samples (noof_samples)

        std::vector<Matrix> sigma = diffusion_batch(X); //compute diffusion matrix(size M, each N x noof_samples)

        // Compute diffusion contribution: N x num_samples
        Matrix diffusion = Matrix::Zero(N, noof_samples);
        for (int n = 0; n < N; ++n) {
            for (int s = 0; s < noof_samples; ++s){
                diffusion(n,s) = (sigma[n].row(s)).dot(dWs[t].col(s));
            }
        }

        X += drift + diffusion;   //update X (matrix N x noof_samples)
    }


    Cumulants result;
    result.mean = mean.colwise().mean();
    result.mean_std = mean_std.colwise().mean();
    result.variance = variance.colwise().mean();
    result.variance_std = variance_std.colwise().mean();

    return result;
}


