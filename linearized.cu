#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <cuda_runtime.h>
#include <thrust/copy.h>
void checkCudaError(cudaError_t err, const char* msg) {
    if (err != cudaSuccess) {
        std::cerr << "CUDA error: " << msg << " - " << cudaGetErrorString(err) << std::endl;
        exit(EXIT_FAILURE);
    }
}

const double alpha = 0.27;
const double C = 1.0;
double h_H = 1.001;
const double phi_start = 0.0;
const double phi_end = M_PI;

const double h = 0.001;
const int steps = (phi_end - phi_start) / h;

#define RK4

//const int steps = 100000;
//const double h = (phi_end - phi_start) / steps;

__device__ void rhs(double varphi, double Q, double Phi, double k, double& dQ, double& dPhi, double H) {
    double denom = H - sin(2.0 * varphi);
    double term1 = (k * k) / (alpha * denom);
    double term2 = 2.0 * cos(2.0 * varphi) / denom;

    dQ = term1 * (C * Phi - alpha * C * Q) + (term2 * Phi / alpha);
    dPhi = term1 * (-alpha * C * Phi - C * Q) - (term2 * Phi);
}



__device__ void euler_solver(double k, double& final_Q, double& final_Phi, double initial_Q, double initial_Phi, double H) {
    double Q = initial_Q, Phi = initial_Phi, varphi = phi_start;
    double k1_Q, k1_Phi;

    for (int i = 0; i < steps; ++i) {
        rhs(varphi, Q, Phi, k, k1_Q, k1_Phi, H);

        Q = Q + h * k1_Q;
        Phi = Phi + h * k1_Phi;
        varphi += h;
    }

    final_Q = Q;
    final_Phi = Phi;
}

__device__ void rk4_solver(double k, double& final_Q, double& final_Phi, double initial_Q, double initial_Phi, double H) {
    double Q = initial_Q, Phi = initial_Phi, varphi = phi_start;
    double k1_Q, k1_Phi, k2_Q, k2_Phi, k3_Q, k3_Phi, k4_Q, k4_Phi;
    double Qt, Phit;

    for (int i = 0; i < steps; ++i) {
        rhs(varphi, Q, Phi, k, k1_Q, k1_Phi, H);

        Qt = Q + 0.5 * h * k1_Q;
        Phit = Phi + 0.5 * h * k1_Phi;
        rhs(varphi + 0.5 * h, Qt, Phit, k, k2_Q, k2_Phi, H);

        Qt = Q + 0.5 * h * k2_Q;
        Phit = Phi + 0.5 * h * k2_Phi;
        rhs(varphi + 0.5 * h, Qt, Phit, k, k3_Q, k3_Phi, H);

        Qt = Q + h * k3_Q;
        Phit = Phi + h * k3_Phi;
        rhs(varphi + h, Qt, Phit, k, k4_Q, k4_Phi, H);

        Q += h / 6.0 * (k1_Q + 2*k2_Q + 2*k3_Q + k4_Q);
        Phi += h / 6.0 * (k1_Phi + 2*k2_Phi + 2*k3_Phi + k4_Phi);
        varphi += h;
    }

    final_Q = Q;
    final_Phi = Phi;
}

__device__ void rk6_solver(double k, double& final_Q, double& final_Phi, double initial_Q, double initial_Phi, double H) {
    double Q = initial_Q, Phi = initial_Phi, varphi = phi_start;

    double a2 = 1.0 / 3.0, a3 = 2.0 / 3.0, a4 = 1.0, a5 = 0.5, a6 = 1.0;
    double b1 = 1.0 / 6.0, b2 = 0.0, b3 = 0.0, b4 = 2.0 / 3.0, b5 = 1.0 / 6.0, b6 = 0.0;

    double k_Q[6], k_Phi[6];
    double Qt, Phit;

    for (int i = 0; i < steps; ++i) {
        rhs(varphi, Q, Phi, k, k_Q[0], k_Phi[0],H);

        Qt = Q + h * a2 * k_Q[0];
        Phit = Phi + h * a2 * k_Phi[0];
        rhs(varphi + h * a2, Qt, Phit, k, k_Q[1], k_Phi[1],H);

        Qt = Q + h * a3 * k_Q[1];
        Phit = Phi + h * a3 * k_Phi[1];
        rhs(varphi + h * a3, Qt, Phit, k, k_Q[2], k_Phi[2],H);

        Qt = Q + h * a4 * k_Q[2];
        Phit = Phi + h * a4 * k_Phi[2];
        rhs(varphi + h * a4, Qt, Phit, k, k_Q[3], k_Phi[3],H);

        Qt = Q + h * a5 * k_Q[3];
        Phit = Phi + h * a5 * k_Phi[3];
        rhs(varphi + h * a5, Qt, Phit, k, k_Q[4], k_Phi[4],H);

        Qt = Q + h * a6 * k_Q[4];
        Phit = Phi + h * a6 * k_Phi[4];
        rhs(varphi + h * a6, Qt, Phit, k, k_Q[5], k_Phi[5],H);

        Q += h * (b1 * k_Q[0] + b4 * k_Q[3] + b5 * k_Q[4]);
        Phi += h * (b1 * k_Phi[0] + b4 * k_Phi[3] + b5 * k_Phi[4],H);

        varphi += h;
    }

    final_Q = Q;
    final_Phi = Phi;
}

__device__ void rk45_solver(
    double k, double& final_Q, double& final_Phi,
    double initial_Q, double initial_Phi, double H
) {
    double Q = initial_Q, Phi = initial_Phi;
    double varphi = phi_start;

    double k_Q[7], k_Phi[7];
    double Qt, Phit;

    for (int i = 0; i < steps; ++i) {
        // c values (nodes)
        double c2 = 1.0 / 5.0;
        double c3 = 3.0 / 10.0;
        double c4 = 4.0 / 5.0;
        double c5 = 8.0 / 9.0;
        double c6 = 1.0;
        double c7 = 1.0;

        // a coefficients
        // Stage 1 (no a's)
        rhs(varphi, Q, Phi, k, k_Q[0], k_Phi[0]);

        // Stage 2
        Qt = Q + h * (1.0 / 5.0) * k_Q[0];
        Phit = Phi + h * (1.0 / 5.0) * k_Phi[0];
        rhs(varphi + c2 * h, Qt, Phit, k, k_Q[1], k_Phi[1]);

        // Stage 3
        Qt = Q + h * (3.0/40.0 * k_Q[0] + 9.0/40.0 * k_Q[1]);
        Phit = Phi + h * (3.0/40.0 * k_Phi[0] + 9.0/40.0 * k_Phi[1]);
        rhs(varphi + c3 * h, Qt, Phit, k, k_Q[2], k_Phi[2]);

        // Stage 4
        Qt = Q + h * (44.0/45.0 * k_Q[0] - 56.0/15.0 * k_Q[1] + 32.0/9.0 * k_Q[2]);
        Phit = Phi + h * (44.0/45.0 * k_Phi[0] - 56.0/15.0 * k_Phi[1] + 32.0/9.0 * k_Phi[2]);
        rhs(varphi + c4 * h, Qt, Phit, k, k_Q[3], k_Phi[3]);

        // Stage 5
        Qt = Q + h * (19372.0/6561.0 * k_Q[0] - 25360.0/2187.0 * k_Q[1]
                   + 64448.0/6561.0 * k_Q[2] - 212.0/729.0 * k_Q[3]);
        Phit = Phi + h * (19372.0/6561.0 * k_Phi[0] - 25360.0/2187.0 * k_Phi[1]
                   + 64448.0/6561.0 * k_Phi[2] - 212.0/729.0 * k_Phi[3]);
        rhs(varphi + c5 * h, Qt, Phit, k, k_Q[4], k_Phi[4]);

        // Stage 6
        Qt = Q + h * (9017.0/3168.0 * k_Q[0] - 355.0/33.0 * k_Q[1]
                   + 46732.0/5247.0 * k_Q[2] + 49.0/176.0 * k_Q[3]
                   - 5103.0/18656.0 * k_Q[4]);
        Phit = Phi + h * (9017.0/3168.0 * k_Phi[0] - 355.0/33.0 * k_Phi[1]
                   + 46732.0/5247.0 * k_Phi[2] + 49.0/176.0 * k_Phi[3]
                   - 5103.0/18656.0 * k_Phi[4]);
        rhs(varphi + c6 * h, Qt, Phit, k, k_Q[5], k_Phi[5]);

        // Stage 7
        Qt = Q + h * (35.0/384.0 * k_Q[0] + 0.0 * k_Q[1] + 500.0/1113.0 * k_Q[2]
                   + 125.0/192.0 * k_Q[3] - 2187.0/6784.0 * k_Q[4]
                   + 11.0/84.0 * k_Q[5]);
        Phit = Phi + h * (35.0/384.0 * k_Phi[0] + 0.0 * k_Phi[1] + 500.0/1113.0 * k_Phi[2]
                   + 125.0/192.0 * k_Phi[3] - 2187.0/6784.0 * k_Phi[4]
                   + 11.0/84.0 * k_Phi[5]);
        rhs(varphi + c7 * h, Qt, Phit, k, k_Q[6], k_Phi[6]);

        // 5th-order solution (you can also compute 4th for adaptive step control)
        Q += h * (35.0/384.0 * k_Q[0] + 500.0/1113.0 * k_Q[2]
                + 125.0/192.0 * k_Q[3] - 2187.0/6784.0 * k_Q[4]
                + 11.0/84.0 * k_Q[5]);

        Phi += h * (35.0/384.0 * k_Phi[0] + 500.0/1113.0 * k_Phi[2]
                + 125.0/192.0 * k_Phi[3] - 2187.0/6784.0 * k_Phi[4]
                + 11.0/84.0 * k_Phi[5]);

        varphi += h;
    }

    final_Q = Q;
    final_Phi = Phi;
}

__device__ void ode_solver(double k, double& final_Q, double& final_Phi, double initial_Q, double initial_Phi, double H)
{
    #ifdef RK4
    rk4_solver(k, final_Q, final_Phi, initial_Q, initial_Phi, H);
    #elif defined(EULER)
    euler_solver(k, final_Q, final_Phi, initial_Q, initial_Phi, H);
    #elif defined(RK6)
    rk6_solver(k, final_Q, final_Phi, initial_Q, initial_Phi, H);
    #elif defined(RK45)
    rk45_solver(k, final_Q, final_Phi, initial_Q, initial_Phi, H);
    #endif
}



__device__ void eigenvalues_magnitudes_2x2(double a, double b, double c, double d, double* mag1, double* mag2) {
    double trace = a + d;
    double det = a * d - b * c;
    double discriminant = trace * trace - 4 * det;

    if (discriminant >= 0) {
        // Real eigenvalues
        double sqrt_disc = sqrt(discriminant);
        double lambda1 = (trace + sqrt_disc) / 2.0;
        double lambda2 = (trace - sqrt_disc) / 2.0;
        *mag1 = fabs(lambda1);
        *mag2 = fabs(lambda2);
    } else {
        // Complex conjugate eigenvalues
        double real_part = trace / 2.0;
        double imag_part = sqrt(-discriminant) / 2.0;
        double magnitude = sqrt(real_part * real_part + imag_part * imag_part);
        *mag1 = magnitude;
        *mag2 = magnitude;
    }
}



__global__ void solve_all_grid(
    const double* __restrict__ k_vals, const double* __restrict__ H_vals, 
    double* Q_out, double* Phi_out, int Nk, int Nh) 
{
    
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int idy = blockIdx.y * blockDim.y + threadIdx.y;
    
    double Q_ini,Phi_ini;
    if (idx < Nk && idy < Nh) 
    {
        double k = k_vals[idx];
        double H = H_vals[idy];

    	double a, b, c, d;
    
    	Q_ini = 1.0; Phi_ini = 0.0;
        ode_solver(k, a, b, Q_ini, Phi_ini, H);
    	
    	Q_ini = 0.0; Phi_ini = 1.0;
        ode_solver(k, c, d, Q_ini, Phi_ini, H);
    
    	double lambda1, lambda2;
    	eigenvalues_magnitudes_2x2(a, b, c, d, &lambda1, &lambda2);
    
        Q_out[idx + Nk*idy] = lambda1;
        Phi_out[idx + Nk*idy] = lambda2;
        //printf("Thread (%d, %d) processing k=%f, H=%f, (%f, %f)\n", idx, idy, k_vals[idx], H_vals[idy], lambda1, lambda2);
    }
}

int main(int argc, char **argv) {

    double hstart, hend;
    double kstart, kend;
    int Nh, Nk;	

    if(argc>1) hstart = atof(argv[1]);
    if(argc>2) hend = atof(argv[2]);
    if(argc>3) Nh = atoi(argv[3]);
    if(argc>4) kstart = atof(argv[4]);
    if(argc>5) kend = atof(argv[5]);
    if(argc>6) Nk = atoi(argv[6]);
    	
    thrust::host_vector<double> h_k_vals(Nk);
    thrust::host_vector<double> h_H_vals(Nh);

    for (int i = 0; i < Nk; ++i)
        h_k_vals[i] = kstart + i * (kend - kstart) / (Nk - 1);

    for (int i = 0; i < Nh; ++i)
        h_H_vals[i] = hstart + i * (hend - hstart) / (Nh - 1);

    thrust::device_vector<double> d_k_vals = h_k_vals;
    thrust::device_vector<double> d_H_vals = h_H_vals;
    thrust::device_vector<double> d_Q_out(Nk*Nh), d_Phi_out(Nk*Nh);

    //assert(h<0.001);

    dim3 blockSize = dim3(16,16);
    int numBlocks_k = (Nk + blockSize.x - 1) / blockSize.x;
    int numBlocks_H = (Nh + blockSize.y - 1) / blockSize.y;
    dim3 numBlocks = dim3(numBlocks_k, numBlocks_H);

    FILE* f = fopen("rk4_k_sweep.txt", "w");
    //fprintf(f, "# k H  |lambda1| |lambda2|\n");

    thrust::host_vector<double> h_Q_out(Nk*Nh);
    thrust::host_vector<double> h_Phi_out(Nk*Nh);

    std::cout << "Nk: " << Nk << ", Nh: " << Nh << std::endl;
    std::cout << "k range: [" << kstart << ", " << kend << "], H range: [" << hstart << ", " << hend << "]" << std::endl;
    std::cout << "Block size: " << blockSize.x << "x" << blockSize.y << ", Number of blocks: " << numBlocks.x << "x" << numBlocks.y << std::endl;
    printf("Starting computation...\n");
    
 // 1. Create events
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    // 2. Record start
    cudaEventRecord(start);
    
    // 3. Launch the kernel
    solve_all_grid<<<numBlocks, blockSize>>>(
                                            thrust::raw_pointer_cast(d_k_vals.data()),
                                            thrust::raw_pointer_cast(d_H_vals.data()),
                                            thrust::raw_pointer_cast(d_Q_out.data()),
                                            thrust::raw_pointer_cast(d_Phi_out.data()), Nk, Nh
                                        );
    checkCudaError(cudaGetLastError(), "Kernel launch failed");
    
    // 4. Record stop
    cudaEventRecord(stop);

    // 5. Wait for the kernel to finish
    cudaEventSynchronize(stop);

    // 6. Compute elapsed time
    float milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, stop);

    std::cout << "Kernel execution time: " << milliseconds << " ms\n";
    
    // 7. Destroy events
    cudaEventDestroy(start);
    cudaEventDestroy(stop);    

    thrust::copy(d_Q_out.begin(),d_Q_out.end(), h_Q_out.begin());
    thrust::copy(d_Phi_out.begin(),d_Phi_out.end(), h_Phi_out.begin());

    for (int j = 0; j < Nh; ++j){
        for (int i = 0; i < Nk; ++i) {
            fprintf(f, "%lf %lf %lf %lf\n", h_k_vals[i], h_H_vals[j], h_Q_out[i + Nk * j], h_Phi_out[i + Nk * j]);
        }
        fprintf(f, "\n");
    }        
 
    fclose(f);
    printf("Results saved to rk4_k_sweep.txt\n");
    return 0;
}
