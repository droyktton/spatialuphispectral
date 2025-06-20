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
const double c = 1.0;
double h_H = 1.001;
const double phi_start = 0.0;
const double phi_end = M_PI;

const double h = 0.0001;
const int steps = (phi_end - phi_start) / h;

//const int steps = 100000;
//const double h = (phi_end - phi_start) / steps;

__device__ void rhs(double varphi, double Q, double Phi, double k, double& dQ, double& dPhi, double H) {
    double denom = H - sin(2.0 * varphi);
    double term1 = (k * k) / (alpha * denom);
    double term2 = 2.0 * cos(2.0 * varphi) / denom;

    dQ = term1 * (c * Phi - alpha * c * Q) + (term2 * Phi / alpha);
    dPhi = term1 * (-alpha * c * Phi - c * Q) - (term2 * Phi);
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
            rk4_solver(k, a, b, Q_ini, Phi_ini, H);
    	
    	Q_ini = 0.0; Phi_ini = 1.0;
            rk4_solver(k, c, d, Q_ini, Phi_ini, H);
    
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
    
    // Launch the kernel
    solve_all_grid<<<numBlocks, blockSize>>>(
                                            thrust::raw_pointer_cast(d_k_vals.data()),
                                            thrust::raw_pointer_cast(d_H_vals.data()),
                                            thrust::raw_pointer_cast(d_Q_out.data()),
                                            thrust::raw_pointer_cast(d_Phi_out.data()), Nk, Nh
                                        );
    checkCudaError(cudaGetLastError(), "Kernel launch failed");

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
