// Libraries
#include "../include/SDE_Solver.hpp"


///////////////////////////   HELPER FUNCTIONS  ////////////////////////////////////

//array<Matrix, noof_timesteps>
// Generate precomputed Brownian increments (shape: time_steps x M x num_samples)
std::vector<Matrix> generateWienerIncrements(const int& M, const int& noof_timesteps, const int& noof_MC_samples, const double& dt) {
    std::mt19937 random_number_generator(std::random_device{}());  //define random number generator object of type mt19937
    std::normal_distribution<> normal_distribution(0.0, 1.0);      //define normal distribution object

    std::vector<Matrix> dWs(noof_timesteps); //generate noof_timesteps matrices with dimension (M x noof_MC_samples)
    for (int t = 0; t < noof_timesteps; ++t) {    //loop over timesteps
        Matrix dW = Matrix::Zero(M, noof_MC_samples);          //dW: contains Wiener increments (M x noof_MC_samples)
        for (int m = 0; m < M; ++m) {             // loop over Wiener dimensions
            for (int s = 0; s < noof_MC_samples; ++s) {   // loop over samples
                dW(m, s) = normal_distribution(random_number_generator);  //calculate Wiener increment
            }
        }
        dWs[t] = dW* std::sqrt(dt);
    }
    return dWs;
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



double integrate_simpson(const Eigen::VectorXd& x, const Eigen::VectorXd& y) {
    // Check input validity
    if (x.size() != y.size()) {
        throw std::invalid_argument("x and y vectors must have the same length.");
    }
    if (x.size() < 2) {
        throw std::invalid_argument("At least 2 points are required.");
    }

    // Check uniform spacing (with floating-point tolerance)
    const double dx = x(1) - x(0);
    const double tolerance = 1e-10;
    for (int i = 1; i < x.size(); ++i) {
        if (std::abs((x(i) - x(i-1)) - dx) > tolerance) {
            throw std::invalid_argument("x must be uniformly spaced.");
        }
    }

    double integral = 0.0;

    // Use Simpson's 1/3 rule for all complete triplets (requires odd number of points)
    int n_simpson = (x.size() % 2 == 1) ? x.size() : x.size() - 1;
    if (n_simpson >= 3) {
        double sum = y(0) + y(n_simpson-1);  // First and last terms in Simpson's range
        for (int i = 1; i < n_simpson - 1; ++i) {
            sum += (i % 2 == 1) ? 4.0 * y(i) : 2.0 * y(i);
        }
        integral += sum * dx / 3.0;
    }

    // Add trapezoidal rule for the last segment if number of points was even
    if (x.size() % 2 == 0 && x.size() >= 3) {
        integral += (y(y.size()-2) + y(y.size()-1)) * dx / 2.0;
    }

    return integral;
}



Vector thomas_algorithm(
    const Vector& sub_diagonal,   // a (size N-1) [L_{i,i-1} for i=2,...,N]
    const Vector& diagonal,       // b (size N)   [L_{i,i} for i=1,...,N]
    const Vector& super_diagonal, // c (size N-1) [L_{i,i+1} for i=1,...,N-1]
    const Vector& rhs             // d (size N)   [right-hand side]
) {
    const int n = diagonal.size();

    // Input validation
    if (sub_diagonal.size() != n - 1 || super_diagonal.size() != n - 1 || rhs.size() != n) {
        throw std::invalid_argument("Invalid input sizes.");
    }

    if (n == 0) return Vector();
    if (n == 1) return Vector::Constant(1, rhs[0] / diagonal[0]);

    Vector c_prime(n - 1);  // Modified super-diagonal
    Vector d_prime(n);      // Modified right-hand side
    Vector x(n);            // Solution

    // Forward elimination
    c_prime[0] = super_diagonal[0] / diagonal[0];
    d_prime[0] = rhs[0] / diagonal[0];

    for (int i = 1; i < n - 1; ++i) {
        double denom = diagonal[i] - sub_diagonal[i - 1] * c_prime[i - 1];
        if (std::abs(denom) < 1e-15) {
            throw std::runtime_error("Zero pivot encountered at index " + std::to_string(i));
        }
        c_prime[i] = super_diagonal[i] / denom;
        d_prime[i] = (rhs[i] - sub_diagonal[i - 1] * d_prime[i - 1]) / denom;
    }

    // Last row (no super-diagonal)
    double denom = diagonal[n - 1] - sub_diagonal[n - 2] * c_prime[n - 2];
    if (std::abs(denom) < 1e-15) {
        throw std::runtime_error("Zero pivot encountered at last row.");
    }
    d_prime[n - 1] = (rhs[n - 1] - sub_diagonal[n - 2] * d_prime[n - 2]) / denom;

    // Back substitution
    x[n - 1] = d_prime[n - 1];
    for (int i = n - 2; i >= 0; --i) {
        x[i] = d_prime[i] - c_prime[i] * x[i + 1];
    }

    return x;
}



///////////////////////////   SDE SOLVERS (LANGEVIN)  ////////////////////////////////////

//general SDE Euler Maruyama Solver that returns mean and variance of an ensemble-averaged trajectory
Cumulants SDE_EulerMaruyama_Solver(
    const int& N, const int& M, //N: state vector dimension, M: Wiener process dimension
    const Vector& X0,   //X0: initial state
    const DriftFunction& drift_batch,   
    const DiffusionFunction& diffusion_batch,
    const double& T, 
    const int& noof_timesteps, const int& noof_MC_samples) {

    Cumulants result;
    result.mean = Matrix::Ones(N, noof_timesteps);  //store mean
    result.variance = Matrix::Ones(N, noof_timesteps); //store variance
    result.mean_std = Matrix::Ones(N, noof_timesteps);  //store std of the mean
    result.variance_std = Matrix::Ones(N, noof_timesteps); //store std of the variance

    Matrix X_centered = Matrix::Ones(N, noof_timesteps);  //store centered mean
    Matrix curtosis = Matrix::Ones(N, noof_timesteps); //store curtosis

    double dt = T / noof_timesteps;

    std::vector<Matrix> dWs = generateWienerIncrements(M, noof_timesteps, noof_MC_samples, dt);  //pre-compute Wiener increments

    Matrix X = X0.replicate(1, noof_MC_samples);  //matrix of initial states

    for (int t = 0; t < noof_timesteps; ++t) {  //loop over noof_timesteps
        
        //Write means first before update(X0)
        result.mean.col(t) = X.rowwise().mean();
        X_centered = X.colwise() - result.mean.col(t);
        //variance.col(t) = (X.colwise() - mean.col(t)).array().square().rowwise().sum()/(noof_MC_samples-1);
        result.variance.col(t) = X_centered.array().square().rowwise().sum() / (noof_MC_samples - 1);

        result.mean_std.col(t) = (result.variance.col(t)/noof_MC_samples).array().sqrt();
        curtosis.col(t) = X_centered.array().pow(4).rowwise().sum() / noof_MC_samples;
        result.variance_std.col(t) = ((curtosis.col(t).array() - result.variance.col(t).array().square()) / noof_MC_samples).sqrt();
        //variance_std.col(t) = ((curtosis.col(t).colwise() - variance.col(t).array().square().matrix())/noof_MC_samples).array().sqrt();

    
        double time = t * dt; //update time

        Matrix drift = drift_batch(X) * dt;  //compute drift for all samples (noof_MC_samples)

        std::vector<Matrix> sigma = diffusion_batch(X); //compute diffusion matrix(size M, each N x noof_MC_samples)

        // Compute diffusion contribution: N x num_samples
        Matrix diffusion = Matrix::Zero(N, noof_MC_samples);
        for (int n = 0; n < N; ++n) {
            for (int s = 0; s < noof_MC_samples; ++s){
                diffusion(n,s) = (sigma[n].row(s)).dot(dWs[t].col(s));
            }
        }

        X += drift + diffusion;   //update X (matrix N x noof_MC_samples)
    }
    return result;
}



// ===== Vectorized Monte Carlo Simulation =====
void SDE_Trajectory_EulerMaruyama_Solver(
    const int& N, const int& M, //N: state vector dimension, M: Wiener process dimension
    const Vector& X0,   //X0: initial state
    const DriftFunction& drift_batch,   
    const DiffusionFunction& diffusion_batch,
    const double& T, 
    const int& noof_timesteps, const int& noof_MC_samples) {

    Cumulants Trajectory_results = SDE_EulerMaruyama_Solver(
    N, M, //N: state vector dimension, M: Wiener process dimension
    X0,   //X0: initial state
    drift_batch,   
    diffusion_batch,
    T, 
    noof_timesteps, noof_MC_samples);

    writeMatrixToFile(Vector::LinSpaced(noof_timesteps, 0, T), "time.txt");
    writeMatrixToFile(Trajectory_results.mean.transpose(), "mean.txt");
    writeMatrixToFile(Trajectory_results.mean_std.transpose(), "mean_std.txt");
    writeMatrixToFile(Trajectory_results.variance.transpose(), "variance.txt");
    writeMatrixToFile(Trajectory_results.variance_std.transpose(), "variance_std.txt");
}


Cumulants SDE_Time_Averaged_EulerMaruyama_Solver(
    const int& N, const int& M, //N: state vector dimension, M: Wiener process dimension
    const Vector& X0,   //X0: initial state
    const DriftFunction& drift_batch,   
    const DiffusionFunction& diffusion_batch,
    const double& T, 
    const int& noof_timesteps, const int& noof_MC_samples) {

    Cumulants Trajectory_results = SDE_EulerMaruyama_Solver(
    N, M, //N: state vector dimension, M: Wiener process dimension
    X0,   //X0: initial state
    drift_batch,   
    diffusion_batch,
    T, 
    noof_timesteps, noof_MC_samples);

    Cumulants Time_averaged_results;
    Time_averaged_results.mean = Trajectory_results.mean.rowwise().mean();
    Time_averaged_results.mean_std = Trajectory_results.mean_std.rowwise().mean();
    Time_averaged_results.variance = Trajectory_results.variance.rowwise().mean();
    Time_averaged_results.variance_std = Trajectory_results.variance_std.rowwise().mean();

    return Time_averaged_results;
}


Long_Time_Mean_Velocity SDE_Long_Time_EulerMaruyama_Solver(
    const int& N, const int& M, //N: state vector dimension, M: Wiener process dimension
    const Vector& X0,   //X0: initial state
    const DriftFunction& drift_batch,   
    const DiffusionFunction& diffusion_batch,
    const double& T, 
    const int& noof_timesteps){

    int noof_MC_samples = 1;

    Long_Time_Mean_Velocity result;
    result.mean_velocity = Matrix::Ones(N, 1);  //store mean
    result.mean_velocity_std = Matrix::Ones(N, 1); //store variance


    double dt = T / noof_timesteps;   //time step size

    int noof_time_samples = 1e3;
    int time_sample_index;
    Matrix time_samples = Matrix::Ones(N, 2*noof_time_samples);  //stores a sample array of X in equidistant time 
    Matrix velocity_samples = Matrix::Ones(N, noof_time_samples);

    std::vector<Matrix> dWs = generateWienerIncrements(M, noof_timesteps, noof_MC_samples, dt);  //pre-compute Wiener increments

    Matrix X = X0.replicate(1, noof_MC_samples);  //matrix of initial states

    for (int t = 0; t < noof_timesteps; ++t) {  //loop over noof_timesteps

        if ((t*2*noof_time_samples)%noof_timesteps == 0){
            time_sample_index = t*2*noof_time_samples/noof_timesteps;
            time_samples.col(time_sample_index) = X;
        }

        double time = t * dt; //update time

        Matrix drift = drift_batch(X) * dt;  //compute drift for all samples (noof_MC_samples)
        std::vector<Matrix> sigma = diffusion_batch(X); //compute diffusion matrix(size M, each N x noof_MC_samples)

        // Compute diffusion contribution: N x num_samples
        Matrix diffusion = Matrix::Zero(N, noof_MC_samples);
        for (int n = 0; n < N; ++n) {
            for (int s = 0; s < noof_MC_samples; ++s){
                diffusion(n,s) = (sigma[n].row(s)).dot(dWs[t].col(s));
            }
        }

        X += drift + diffusion;   //update X (matrix N x noof_MC_samples)
    }
    for (int i = 0; i < noof_time_samples; ++i){
        velocity_samples.col(i) = (time_samples.col(i+ noof_time_samples).array() - time_samples.col(i).array())*2.0/T;
    }
    result.mean_velocity = velocity_samples.rowwise().mean();
    Matrix time_sample_differences_centered = velocity_samples- result.mean_velocity.replicate(1, noof_time_samples);
    result.mean_velocity_std = (((time_sample_differences_centered.array().square()).rowwise().sum())/(noof_time_samples-1)).array().sqrt();
    result.mean_velocity_std = result.mean_velocity_std.array()/std::sqrt(noof_time_samples);
    return result;
}




SDE_1D_result SDE_1D_EulerMaruyama_Solver(
    Fokker_Planck_Drift_Function drift_function, 
    const double& D,
    const double& x0,   //X0: initial state 
    const int& noof_MC_samples,       //Number of samples 
    const int& noof_timesteps,
    const int& update_frequency,
    double dt, 
    double tolerance_CLT){

    SDE_1D_result result;
    vector <double> result_time;
    vector <double> result_mean;
    vector <double> result_variance;
    vector <double> result_mean_velocity;
    vector <double> result_variance_velocity;
    vector <double> result_mean_std;
    vector <double> result_variance_std;
    vector <double> result_mean_velocity_std;
    vector <double> result_variance_velocity_std;
    Vector X_centered = Vector::Zero(noof_MC_samples);
    double fourth_central_moment;
    double variance_uncertainty_squared;
    Vector X = Vector::Ones(noof_MC_samples)*x0;

    std::mt19937 random_number_generator(std::random_device{}());  //define random number generator object of type mt19937
    std::normal_distribution<> normal_distribution(0.0, 1.0);      //define normal distribution object
    double random_number;

    double drift;
    double diffusion;

    double relative_change_mean_velocity;
    double relative_change_variance_velocity;


    for (int n = 0; n < noof_timesteps; ++n){
        for (int i = 0; i < noof_MC_samples; ++i){
            random_number = normal_distribution(random_number_generator);
            drift = drift_function(X(i))*dt;
            diffusion = std::sqrt(2*D*dt)*random_number;
            X(i) +=  drift + diffusion;
        }

        // Store results
        if ((n+1) % update_frequency == 0) {
            int idx = (n+1) / update_frequency;

            

            result_time.push_back((n+1) * dt);

            result_mean.push_back(X.mean());

            X_centered = X.array() - result_mean.back();
            result_variance.push_back(X_centered.array().square().sum()/(noof_MC_samples-1));
            
            result_mean_velocity.push_back(result_mean.back()/result_time.back());
            // std::cout << "mean velocity: " << result_mean_velocity.back() << std::endl;
            result_variance_velocity.push_back(result_variance.back()/result_time.back());

            result_mean_std.push_back(std::sqrt(result_variance.back()/noof_MC_samples));

            fourth_central_moment = X_centered.array().pow(4).sum()/(noof_MC_samples-1);
            variance_uncertainty_squared = (fourth_central_moment - std::pow(result_variance.back(), 2))/(noof_MC_samples-1);
            result_variance_std.push_back(std::sqrt(variance_uncertainty_squared));

            result_mean_velocity_std.push_back(result_mean_std.back()/result_time.back());
            result_variance_velocity_std.push_back(result_variance_std.back()/result_time.back());
            
            std::cout << "mean_velo: " << result_mean_velocity.back() << " variance_velo: " << result_variance_velocity.back() <<std::endl;
            if (idx > 0){
                    double relative_change_mean_velocity = std::abs((result_mean_velocity.back() - result_mean_velocity[result_mean_velocity.size()-2])/result_mean_velocity.back());
                    double relative_change_variance_velocity = std::abs((result_variance_velocity.back() - result_variance_velocity[result_mean_velocity.size()-2])/result_variance_velocity.back());
                    if ((relative_change_mean_velocity < tolerance_CLT) && (relative_change_variance_velocity < tolerance_CLT)){
                        std::cout << "CLT tolerance achieved! mean_velo: " << relative_change_mean_velocity <<  " variance_velo: "  <<relative_change_variance_velocity << std::endl;
                        break;
                    }
            }
        }
    }
    result.time = Eigen::Map<Vector>(result_time.data(), result_time.size());
    result.mean = Eigen::Map<Vector>(result_mean.data(), result_mean.size());
    result.variance = Eigen::Map<Vector>(result_variance.data(), result_variance.size());
    result.mean_velocity = Eigen::Map<Vector>(result_mean_velocity.data(), result_mean_velocity.size());
    result.variance_velocity = Eigen::Map<Vector>(result_variance_velocity.data(), result_variance_velocity.size());
    
    result.mean_std = Eigen::Map<Vector>(result_mean_std.data(), result_mean_std.size());
    result.variance_std = Eigen::Map<Vector>(result_variance_std.data(), result_variance_std.size());
    result.mean_velocity_std = Eigen::Map<Vector>(result_mean_velocity_std.data(), result_mean_velocity_std.size());
    result.variance_velocity_std = Eigen::Map<Vector>(result_variance_velocity_std.data(), result_variance_velocity_std.size());
    return result;
}




///////////////////////////   FOKKER PLANCK FRAMEWORK SOLVERS  ////////////////////////////////////

//Fokker Planck 1D Euler Forward solver 
Fokker_Planck_1D_result Fokker_Planck_FVM_1D_Explicit_Euler_Solver(
    Fokker_Planck_Drift_Function drift_function,
    double D,
    const Vector& x,
    const Vector& p0,
    const int& noof_timesteps,
    const int& update_frequency,
    double dt = -1.0) {

    double dx = x(1) - x(0);
    int N = p0.size();
    
    // Face-centered x values
    Vector x_faces = (x.head(N-1) + x.tail(N-1)) / 2.0;
    
    // Evaluate drift at faces
    Vector a_faces(x_faces.size());
    Vector a_nodes(x.size());
    
    // Initialize solution
    int noof_output_steps = noof_timesteps / update_frequency + 1;
    Fokker_Planck_1D_result result;
    result.time = Vector(noof_output_steps);
    result.mean = Vector(noof_output_steps);
    result.variance = Vector(noof_output_steps);
    result.mean_velocity = Vector(noof_output_steps);
    result.variance_velocity = Vector(noof_output_steps);

    Vector p_current = p0;
    Vector x_centered = x.array() - result.mean(0);


    double tolerance_boundary = 1e-23;  // 10 sigma
    
    // Calculate stable time step
    if (dt <= 0.0) {
        double max_a = a_faces.cwiseAbs().maxCoeff();
        dt = std::min(0.5 * dx / max_a, 0.5 * dx * dx / D) / 10.0;
    }
    
    std::cout << "Time step size: " << dt << std::endl;
    
    // Pre-allocate flux vectors
    Vector j_plus(N-1), j_minus(N-1);
    Vector p_upwind(N-1), p_central(N-1), diffusion_term(N-1);
    
    // Time stepping loop
    for (int n = 0; n < noof_timesteps; ++n) {

        for (int i = 0; i < x_faces.size(); ++i) {
            a_faces(i) = drift_function(x_faces(i));
        }

        // Upwind scheme
        p_upwind = (a_faces.array() > 0).select(
            p_current.head(N-1), 
            p_current.tail(N-1)
        );
        
        p_central = (p_current.tail(N-1) + p_current.head(N-1))/2.0;

        // Diffusion term
        diffusion_term = (p_current.tail(N-1) - p_current.head(N-1)) * D / dx;
        
        // Flux calculation
        // j_plus = a_faces.cwiseProduct(p_upwind) - diffusion_term;
        j_plus = a_faces.cwiseProduct(p_central) - diffusion_term;

        
        // Shift fluxes
        j_minus.setZero();
        j_minus.tail(N-2) = j_plus.head(N-2);
        
        // Boundary conditions
        j_minus(0) = 0.0;
        j_plus(N-2) = 0.0;
        
        // Update probability
        p_current.segment(1, N-2) -= (dt/dx) * (j_plus.tail(N-2) - j_minus.tail(N-2));
        p_current(0) -= (dt/dx) * j_plus(0);
        p_current(N-1) += (dt/dx) * j_minus(N-2);
        
        // Non-negativity and normalization
        p_current = p_current.cwiseMax(0.0);
        // double sum = p_current.sum() * dx;
        double sum = integrate_simpson(x, p_current);
        if (sum > 0) p_current /= sum;

        if ((p_current(N-1) > tolerance_boundary) || (p_current(0) > tolerance_boundary)) {
            std::cout << "Early termination: zero-flux not satisfied" << "\n";
            return result;
        }
        
        // Store results
        if ((n+1) % update_frequency == 0) {
            int idx = (n+1) / update_frequency;
            if (idx < noof_output_steps) {

                result.time(idx) = (n+1) * dt;

                // SIMPSON RULE
                // result.mean(idx) = integrate_simpson(x, p_current.cwiseProduct(x));
                // x_centered = x.array() - result.mean(idx);
                // result.variance(idx) = integrate_simpson(x, p_current.cwiseProduct(x_centered.array().square().matrix()));
                // result.mean_velocity(idx) = integrate_simpson(x_faces, j_plus);
                // result.variance_velocity(idx) = 2*integrate_simpson(x_faces, j_plus.cwiseProduct((x_faces.array()- result.mean(idx)).matrix()));

                result.mean(idx) = integrate_simpson(x, p_current.cwiseProduct(x));
                x_centered = x.array() - result.mean(idx);
                result.variance(idx) = integrate_simpson(x, p_current.cwiseProduct(x_centered.array().square().matrix()));

                for (int i = 0; i < x.size(); ++i) {
                    a_nodes(i) = drift_function(x(i));
                }
                result.mean_velocity(idx) = integrate_simpson(x, p_current.cwiseProduct(a_nodes));
                result.variance_velocity(idx) = 2*integrate_simpson(x, p_current.cwiseProduct(x.cwiseProduct(a_nodes))) + 2*D - 2*result.mean(idx)*result.mean_velocity(idx);

                // std::cout << "Boundary values: [left=" << p_current[0] << ", right=" << p_current[N-1] << "\n";
            }
        }
    }
    
    return result;
}




Fokker_Planck_1D_result Fokker_Planck_FVM_1D_Crank_Nicolson_Solver(
    Fokker_Planck_Drift_Function drift_function,
    double D,
    const Vector& x,
    const Vector& p0,
    const int& noof_timesteps,
    const int& update_frequency,
    double dt = -1.0, 
    double tolerance_CLT = 1e-6) {


    double tolerance_boundary = 1e-23;  // 10 sigma (for checking if the PDF hits a boundary)
    double dx = x(1) - x(0);
    int N = p0.size();
    Vector x_faces = (x.head(N - 1) + x.tail(N - 1)) / 2.0;
    Vector p_current = p0;
    Vector rhs = Vector::Zero(N);

    // Evaluate drift at faces and nodes
    Vector a_faces(x_faces.size());
    Vector a_nodes(x.size());

    for (int i = 0; i < x_faces.size(); ++i) {
            a_faces(i) = drift_function(x_faces(i));
    }
    for (int i = 0; i < x.size(); ++i) {
        a_nodes(i) = drift_function(x(i));
    }

    // Stable time step estimate (if not provided)
    if (dt <= 0.0) {
        double max_a = a_faces.cwiseAbs().maxCoeff();
        dt = std::min(dx / max_a, dx * dx / (2.0*D));
    }

    std::cout << "Crank-Nicolson dt: " << dt << std::endl;


    // Output setup
    int noof_output_steps = noof_timesteps / update_frequency + 1;
    Fokker_Planck_1D_result result;

    vector <double> result_time;
    vector <double> result_mean;
    vector <double> result_variance;
    vector <double> result_mean_velocity;
    vector <double> result_variance_velocity;
    Vector x_centered = Vector::Zero(N);
    Vector dp_dt = Vector::Zero(N);

    
    const double theta = 0.5;  // Crank-Nicolson parameter
    Vector L_diagonal(N), L_sub_diagonal(N-1), L_super_diagonal(N-1);


    L_sub_diagonal = -1.0/dx*(a_faces.array()/2.0 + D/dx);

    L_diagonal(0) = 1.0/dx*(a_faces(1)/2.0 + D/dx);
    L_diagonal.segment(1,N-2) = 1.0/dx * ((a_faces.tail(N-2) - a_faces.head(N-2)).array()/2.0 + 2.0*D/dx); 
    L_diagonal(N-1) = - 1.0/dx*(a_faces(N-2)/2.0 - D/dx);

    L_super_diagonal(0) = 1.0/dx*(a_faces(1)/2.0 - D/dx);
    L_super_diagonal.tail(N-2) = 1.0/dx*(a_faces.tail(N-2).array()/2.0 - D/dx);


    Vector LHS_diagonal(N), LHS_sub_diagonal(N-1), LHS_super_diagonal(N-1);

    LHS_sub_diagonal = dt*theta*L_sub_diagonal;
    LHS_diagonal = Vector::Ones(N) + dt*theta*L_diagonal;
    LHS_super_diagonal = dt*theta*L_super_diagonal;


    for (int n = 0; n < noof_timesteps; ++n) {

        rhs = p_current;
        rhs -= dt*(1.0 - theta)*L_diagonal.cwiseProduct(p_current);
        rhs.tail(N-1) -= dt*(1.0 - theta)*L_sub_diagonal.cwiseProduct(p_current.head(N-1));
        rhs.head(N-1) -= dt*(1.0 - theta)*L_super_diagonal.cwiseProduct(p_current.tail(N-1));

        p_current = thomas_algorithm(LHS_sub_diagonal, LHS_diagonal, LHS_super_diagonal, rhs);

        // Non-negativity and normalization
        p_current = p_current.cwiseMax(0.0);
        // double sum = p_current.sum() * dx;
        double sum = integrate_simpson(x, p_current);
        if (sum > 0) p_current /= sum;

        
        
        // Store results
        if ((n+1) % update_frequency == 0) {
            int idx = (n+1) / update_frequency;
            if (idx < noof_output_steps) {

                if ((p_current(N-1) > tolerance_boundary) || (p_current(0) > tolerance_boundary)) {
                    std::cout << "Early termination: zero-flux at boundary not satisfied" << "\n";
                    return result;
                }

                result_time.push_back((n+1) * dt);

                result_mean.push_back(integrate_simpson(x, p_current.cwiseProduct(x)));
                x_centered = x.array() - result_mean.back();
                result_variance.push_back(integrate_simpson(x, p_current.cwiseProduct(x_centered.array().square().matrix())));
                
                dp_dt = -L_diagonal.cwiseProduct(p_current);
                dp_dt.tail(N-1) -= L_sub_diagonal.cwiseProduct(p_current.head(N-1));
                dp_dt.head(N-1) -= L_super_diagonal.cwiseProduct(p_current.tail(N-1));

                result_mean_velocity.push_back(integrate_simpson(x, x.cwiseProduct(dp_dt)));
                result_variance_velocity.push_back(integrate_simpson(x, x_centered.array().square().matrix().cwiseProduct(dp_dt)));

                if (idx > 0){
                    double relative_change_mean_velocity = std::abs((result_mean_velocity.back() - result_mean_velocity[result_mean_velocity.size()-2])/result_mean_velocity.back());
                    double relative_change_variance_velocity = std::abs((result_variance_velocity.back() - result_variance_velocity[result_mean_velocity.size()-2])/result_variance_velocity.back());
                    if ((relative_change_mean_velocity < tolerance_CLT) && (relative_change_variance_velocity < tolerance_CLT)){
                        std::cout << "CLT tolerance achieved! mean_velo: " << relative_change_mean_velocity <<  " variance_velo: "  <<relative_change_variance_velocity << std::endl;
                        break;
                    }
                }
            }
        }
    }
    result.time = Eigen::Map<Vector>(result_time.data(), result_time.size());
    result.mean = Eigen::Map<Vector>(result_mean.data(), result_mean.size());
    result.variance = Eigen::Map<Vector>(result_variance.data(), result_variance.size());
    result.mean_velocity = Eigen::Map<Vector>(result_mean_velocity.data(), result_mean_velocity.size());
    result.variance_velocity = Eigen::Map<Vector>(result_variance_velocity.data(), result_variance_velocity.size());
    
    return result;
}