#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <cmath>
#include <cstdio>

const double alpha = 1.0;
const double c = 1.0;
__device__ double H;
double h_H = 1.001;
const double phi_start = 0.0;
const double phi_end = 10*M_PI;

const double h = 0.0001;
const int steps = (phi_end - phi_start) / h;

//const int steps = 100000;
//const double h = (phi_end - phi_start) / steps;

__device__ void rhs(double varphi, double Q, double Phi, double k, double& dQ, double& dPhi) {
    double denom = H - sin(2.0 * varphi);
    double term1 = (k * k) / (alpha * denom);
    double term2 = 2.0 * cos(2.0 * varphi) / denom;

    dQ = term1 * (c * Phi - alpha * c * Q) + (term2 * Phi / alpha);
    dPhi = term1 * (-alpha * c * Phi - c * Q) - (term2 * Phi);
}

__device__ void rk4_solver(double k, double& final_Q, double& final_Phi) {
    double Q = 1.0, Phi = 0.0, varphi = phi_start;
    double k1_Q, k1_Phi, k2_Q, k2_Phi, k3_Q, k3_Phi, k4_Q, k4_Phi;
    double Qt, Phit;

    for (int i = 0; i < steps; ++i) {
        rhs(varphi, Q, Phi, k, k1_Q, k1_Phi);

        Qt = Q + 0.5 * h * k1_Q;
        Phit = Phi + 0.5 * h * k1_Phi;
        rhs(varphi + 0.5 * h, Qt, Phit, k, k2_Q, k2_Phi);

        Qt = Q + 0.5 * h * k2_Q;
        Phit = Phi + 0.5 * h * k2_Phi;
        rhs(varphi + 0.5 * h, Qt, Phit, k, k3_Q, k3_Phi);

        Qt = Q + h * k3_Q;
        Phit = Phi + h * k3_Phi;
        rhs(varphi + h, Qt, Phit, k, k4_Q, k4_Phi);

        Q += h / 6.0 * (k1_Q + 2*k2_Q + 2*k3_Q + k4_Q);
        Phi += h / 6.0 * (k1_Phi + 2*k2_Phi + 2*k3_Phi + k4_Phi);
        varphi += h;
    }

    final_Q = Q;
    final_Phi = Phi;
}

__global__ void solve_all_k(const double* __restrict__ k_vals, double* Q_out, double* Phi_out, int N) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < N) {
        double k = k_vals[idx];
        double Q, Phi;
        rk4_solver(k, Q, Phi);
        Q_out[idx] = Q;
        Phi_out[idx] = Phi;
    }
}

int main(int argc, char **argv) {

    if(argc>1) h_H = atof(argv[1]);
    cudaMemcpyToSymbol(H, &h_H, sizeof(double));

    std::cout << "H = " << h_H << std::endl;

    const int N = 1024;
    thrust::host_vector<double> h_k_vals(N);
    double kmin = 0.01, kmax = M_PI;

    for (int i = 0; i < N; ++i)
        h_k_vals[i] = kmin + i * (kmax - kmin) / (N - 1);

    thrust::device_vector<double> d_k_vals = h_k_vals;
    thrust::device_vector<double> d_Q_out(N), d_Phi_out(N);

    //assert(h<0.001);

    int blockSize = 128;
    int numBlocks = (N + blockSize - 1) / blockSize;

    FILE* f = fopen("rk4_k_sweep.txt", "w");
    fprintf(f, "# k Q Phi\n");

    thrust::host_vector<double> h_Q_out(N);
    thrust::host_vector<double> h_Phi_out(N);

    // loop over periods
    for(int n=0;n<5;n++){
      solve_all_k<<<numBlocks, blockSize>>>(thrust::raw_pointer_cast(d_k_vals.data()),
                                            thrust::raw_pointer_cast(d_Q_out.data()),
                                            thrust::raw_pointer_cast(d_Phi_out.data()), N);

      cudaDeviceSynchronize();

      thrust::copy(d_Q_out.begin(),d_Q_out.end(), h_Q_out.begin());
      thrust::copy(d_Phi_out.begin(),d_Phi_out.end(), h_Phi_out.begin());

      for (int i = 0; i < N; ++i)
          fprintf(f, "%.6f %.6f %.6f %d\n", h_k_vals[i], h_Q_out[i], h_Phi_out[i], n);
      //fprintf(f, "\n");
    }

    fclose(f);
    printf("Results saved to rk4_k_sweep.txt\n");
    return 0;
}
