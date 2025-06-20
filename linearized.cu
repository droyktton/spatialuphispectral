#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <cmath>
#include <cstdio>

const double alpha = 0.27;
const double c = 1.0;
__device__ double H;
double h_H = 1.001;
const double phi_start = 0.0;
const double phi_end = M_PI;

const double h = 0.00001;
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



__global__ void solve_all_k(const double* __restrict__ k_vals, double* Q_out, double* Phi_out, int N, double H) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    double Q_ini,Phi_ini;
    if (idx < N) {
        double k = k_vals[idx];
        //double Q, Phi;

	double a, b, c, d;

	Q_ini = 1.0; Phi_ini = 0.0;
        rk4_solver(k, a, b, Q_ini, Phi_ini, H);
	
	Q_ini = 0.0; Phi_ini = 1.0;
        rk4_solver(k, c, d, Q_ini, Phi_ini, H);

	double lambda1, lambda2;
	eigenvalues_magnitudes_2x2(a, b, c, d, &lambda1, &lambda2);

        Q_out[idx] = lambda1;
        Phi_out[idx] = lambda2;
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
    	
    //cudaMemcpyToSymbol(H, &h_H, sizeof(double));
    //std::cout << "H = " << h_H << std::endl;

    //const int N = 1024;
    thrust::host_vector<double> h_k_vals(Nk);
    double kmin = kstart, kmax = kend;

    for (int i = 0; i < Nk; ++i)
        h_k_vals[i] = kmin + i * (kmax - kmin) / (Nk - 1);

    thrust::device_vector<double> d_k_vals = h_k_vals;
    thrust::device_vector<double> d_Q_out(Nk), d_Phi_out(Nk);

    //assert(h<0.001);

    int blockSize = 128;
    int numBlocks = (Nk + blockSize - 1) / blockSize;

    FILE* f = fopen("rk4_k_sweep.txt", "w");
    //fprintf(f, "# k H  |lambda1| |lambda2|\n");

    thrust::host_vector<double> h_Q_out(Nk);
    thrust::host_vector<double> h_Phi_out(Nk);

    // loop over periods
    for(int j=0 ; j<Nh ; ++j){
      
      h_H = hstart + j * (hend-hstart) / (Nh-1);		
      //cudaMemcpyToSymbol(H, &h_H, sizeof(double));

      solve_all_k<<<numBlocks, blockSize>>>(thrust::raw_pointer_cast(d_k_vals.data()),
                                            thrust::raw_pointer_cast(d_Q_out.data()),
                                            thrust::raw_pointer_cast(d_Phi_out.data()), Nk, h_H);

      //cudaDeviceSynchronize();

      thrust::copy(d_Q_out.begin(),d_Q_out.end(), h_Q_out.begin());
      thrust::copy(d_Phi_out.begin(),d_Phi_out.end(), h_Phi_out.begin());

      for (int i = 0; i < Nk; ++i)
          fprintf(f, "%.6f %.6f %.6f %.6f\n", h_k_vals[i], h_H, h_Q_out[i], h_Phi_out[i]);
      fprintf(f, "\n");
    }

    fclose(f);
    printf("Results saved to rk4_k_sweep.txt\n");
    return 0;
}
