#include <thrust/complex.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/reduce.h>
#include <thrust/transform_reduce.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/tuple.h>
#include <thrust/functional.h>

//#include <cufftXt.h>
#include <cufft.h>
#include <cuda_runtime.h>
#include <cmath>
#include <iostream>
#include <fstream>

#ifdef DOUBLE_PRECISION
typedef double REAL;
#else
typedef float REAL;
#endif

using complex = thrust::complex<REAL>;


/*#ifndef N
#define N 32768
#endif
*/
//const int N = 32768;

int h_N;
__constant__ int N;

REAL h_Ba;
__constant__ REAL B_a;

REAL L, dx, dt;
int steps;
complex alpha;
REAL K, N_n;

/*const REAL L = h_N*1.0f;
const REAL dx = L / h_N;
const REAL dt = 0.2f;
const int steps = 500000;

const complex alpha(0.27f, 0.0f);
const REAL K = 0.796f;
const REAL N_n = 0.016f;
*/

#ifndef EPSILON
#define EPSILON 0.0001f
#endif

#define NBINS 100


complex one_particle_solution(REAL h){
    REAL a = alpha.real();
    REAL hw = a*N_n/2.0f; 
    REAL vphi = (h > hw)?(sqrtf((h/hw)*(h/hw)-1.0f)*hw/(1+a*a)):0.0f;  
    REAL vu = (h > hw)?(h/a-sqrtf((h/hw)*(h/hw)-1.0f)*hw/(a+a*a*a)):(h/a);
    complex z(vu, vphi);
    return z;  
}


__global__ void histogramKernel(const float* data, int* bins, int N, int Nbins, float xmin, float xmax) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < N) {
        float x = data[idx];
        int bin = int((x - xmin) / (xmax - xmin) * Nbins);
        if (bin >= 0 && bin < Nbins) {
            atomicAdd(&bins[bin], 1);
        }
    }
}

__global__ void init_wave_numbers(complex* L_k, REAL K, REAL L) {
    int i = threadIdx.x + blockDim.x * blockIdx.x;
    if (i < N) {
        int k = (i <= N/2) ? i : i - N;
        REAL kx = 2 * M_PI * k / L;
        L_k[i] = complex(-K * kx * kx, 0.0f);
    }
}

__global__ void nonlinear_term_kernel(const complex* z, complex* nonlinear, REAL N_n) {
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i < N) {
        REAL phi = -z[i].imag();
        //nonlinear[i] = -complex(0.0f, 1.0f) * (N_n / 2.0f) * sinf(2 * phi); // different in the notes...
        nonlinear[i] = complex(0.0f, 1.0f) * (N_n / 2.0f) * sinf(2 * phi);
    }
}

__global__ void crank_nicholson_update(complex* z_hat, const complex* N_hat,
                                       const complex* L_k,
                                       complex alpha, REAL dt, REAL B_a) {
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i < N) {
        complex i_unit(0,1);
        complex a_plus = alpha + i_unit + 0.5f * dt * L_k[i];
        complex a_minus = alpha + i_unit - 0.5f * dt * L_k[i];
        complex rhs = a_plus * z_hat[i] + dt * N_hat[i];
        if (i == 0) rhs -= dt * B_a * static_cast<REAL>(N);
        z_hat[i] = rhs / a_minus;
    }
}

__global__ void normalize(complex* z) {
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i < N) {
        z[i] /= static_cast<REAL>(N);
    }
}

thrust::tuple<complex,complex> roughness(thrust::device_vector<complex> &z)
{
    size_t N_ = z.size();
    complex zcm = thrust::reduce(z.begin(), z.end());
    zcm *= 1.0f/N_;

    complex zcm2 = thrust::transform_reduce(
        z.begin(),
        z.end(),
        [zcm]__device__ __host__ (complex z) {
          REAL realdiff=z.real()-zcm.real();
          REAL imagdiff=z.imag()-zcm.imag();
          return complex(
            realdiff*realdiff,
            imagdiff*imagdiff
          );
        },
        complex(0.0f, 0.0f),
        thrust::plus<complex>());

    zcm2 *= 1.0f/N_;

    return make_tuple(zcm2, zcm);
}

class Cuerda
{
    public:
    Cuerda(){
      z.resize(h_N);
      z_hat.resize(h_N);
      nonlinear.resize(h_N);
      L_k.resize(h_N);
      z_prev.resize(h_N);
      dzdt.resize(h_N);
      //zaux.resize(N);
      histogram_udot.resize(NBINS);
      histogram_phidot.resize(NBINS);

      init();

      #ifdef DEBUG
      std::cout << "N=" << h_N << std::endl;
      #endif
    };
    
    REAL distance_conf_fourier(Cuerda &c, int nc){
        complex zcm1=complex(0.f,0.f);
        complex zcm2=complex(0.f,0.f);
        
        REAL dist =
        thrust::transform_reduce
        (
          thrust::make_zip_iterator(thrust::make_tuple(z_hat.begin()+1,c.z_hat.begin()+1)),
          thrust::make_zip_iterator(thrust::make_tuple(z_hat.begin()+nc+1,c.z_hat.begin()+nc+1)),
          [zcm1,zcm2]__device__ __host__ (thrust::tuple<complex,complex> tuplez)
          {
            complex z1 = thrust::get<0>(tuplez)-zcm1;
            complex z2 = thrust::get<1>(tuplez)-zcm1;
            return powf(z1.real() - z2.real(),2.0)+powf(z1.imag() - z2.imag(),2.0);
            //return powf(z1.real() - z2.real(),2.0);
          },
          REAL(0.0f),
          thrust::plus<REAL>()
        );

        return dist/h_N;
    }

    REAL distance_conf(Cuerda &c){
      
      complex zcm1=complex(0.f,0.f);
      complex zcm2=complex(0.f,0.f);

      /*complex zcm1 = thrust::reduce(z.begin(), z.end());
      complex zcm2 = thrust::reduce(c.z.begin(), c.z.end());
      zcm1 *= 0.0f/N;
      zcm2 *= 0.0f/N;*/


      REAL dist =
      thrust::transform_reduce
      (
        thrust::make_zip_iterator(thrust::make_tuple(z.begin(),c.z.begin())),
        thrust::make_zip_iterator(thrust::make_tuple(z.end(),c.z.end())),
        [zcm1,zcm2]__device__ __host__ (thrust::tuple<complex,complex> tuplez)
        {
          complex z1 = thrust::get<0>(tuplez);
          complex z2 = thrust::get<1>(tuplez);
          return powf(z1.real() - z2.real(),2.0)+powf(z1.imag() - z2.imag(),2.0);
          //return powf(z1.real() - z2.real(),2.0);
        },
        REAL(0.0f),
        thrust::plus<REAL>()
      );

      return dist/h_N;
    }

    void copy_conf(Cuerda &c)
    {
      thrust::copy(c.z.begin(),c.z.end(),z.begin());
      thrust::copy(c.z_hat.begin(),c.z_hat.end(),z_hat.begin());
      thrust::copy(c.nonlinear.begin(),c.nonlinear.end(),nonlinear.begin());
      thrust::copy(c.L_k.begin(),c.L_k.end(),L_k.begin());
    }

    void perturb_conf(REAL epsilon)
    {
      thrust::host_vector<complex> hz(h_N);
      complex avhz;
      for (int i = 0; i < h_N; ++i) {
        hz[i] = epsilon*complex(rand()*1.0f/RAND_MAX, rand()*1.0f/RAND_MAX);
        avhz += hz[i];
      }
      avhz *= 1.0f/h_N;
      for(int i=0; i<h_N; i++) z[i] += (hz[i]-avhz);
    }

    ~Cuerda(){
      cufftDestroy(plan);
    };

    void init(){
      #ifdef DOUBLE_PRECISION
      cufftPlan1d(&plan, h_N, CUFFT_Z2Z, 1);
      #else 
      cufftPlan1d(&plan, h_N, CUFFT_C2C, 1);
      #endif
      
      srand(42);
      // Initial condition: z = cos(x)
      //thrust::host_vector<complex> z0(h_N);
      
      /*for (int i = 0; i < h_N; ++i) {
          REAL x = i * dx;
          //z0[i] = complex(0.0f, 0.0f);
          z0[i] = 0.0001*complex(rand()*1.0f/RAND_MAX, rand()*1.0f/RAND_MAX);
          //if(i<10) std::cout << z0[i] << std::endl;
      }
      z = z0;*/
      
      thrust::fill(z.begin(), z.end(), complex(0.0f, 0.0f));
      perturb_conf(EPSILON);
      thrust::copy(z.begin(), z.end(), z_prev.begin());
      thrust::fill(dzdt.begin(), dzdt.end(), complex(0.0f, 0.0f));
      
      // Init linear operator
      init_wave_numbers<<<(h_N+255)/256, 256>>>(thrust::raw_pointer_cast(L_k.data()), K, L);

      #ifdef DEBUG
      std::cout << "K=" << K << std::endl;
      std::cout << "L=" << L << std::endl;
      std::cout << "dx=" << dx << std::endl;
      std::cout << "dt=" << dt << std::endl;
      std::cout << "steps=" << steps << std::endl;
      std::cout << "alpha=" << alpha << std::endl;
      std::cout << "N_n=" << N_n << std::endl;
      std::cout << "B_a=" << h_Ba << std::endl;
      std::cout << "N=" << h_N << std::endl;
      #endif
    }

    void step(){
       nonlinear_term_kernel<<<(h_N+255)/256, 256>>>(
            thrust::raw_pointer_cast(z.data()),
            thrust::raw_pointer_cast(nonlinear.data()), N_n);

        #ifdef DOUBLE_PRECISION
        cufftExecZ2Z(plan,
            reinterpret_cast<cufftDoubleComplex*>(thrust::raw_pointer_cast(z.data())),
            reinterpret_cast<cufftDoubleComplex*>(thrust::raw_pointer_cast(z_hat.data())),
            CUFFT_FORWARD);
        #else
        cufftExecC2C(plan,
            reinterpret_cast<cufftComplex*>(thrust::raw_pointer_cast(z.data())),
            reinterpret_cast<cufftComplex*>(thrust::raw_pointer_cast(z_hat.data())),
            CUFFT_FORWARD);
        #endif


        #ifdef DOUBLE_PRECISION
        cufftExecZ2Z(plan,
            reinterpret_cast<cufftDoubleComplex*>(thrust::raw_pointer_cast(nonlinear.data())),
            reinterpret_cast<cufftDoubleComplex*>(thrust::raw_pointer_cast(nonlinear.data())),
            CUFFT_FORWARD);
        #else
        cufftExecC2C(plan,
            reinterpret_cast<cufftComplex*>(thrust::raw_pointer_cast(nonlinear.data())),
            reinterpret_cast<cufftComplex*>(thrust::raw_pointer_cast(nonlinear.data())),
            CUFFT_FORWARD);
        #endif
      
        crank_nicholson_update<<<(h_N+255)/256, 256>>>(
            thrust::raw_pointer_cast(z_hat.data()),
            thrust::raw_pointer_cast(nonlinear.data()),
            thrust::raw_pointer_cast(L_k.data()),
            alpha, dt, h_Ba);

        #ifdef DOUBLE_PRECISION
        cufftExecZ2Z(plan,
            reinterpret_cast<cufftDoubleComplex*>(thrust::raw_pointer_cast(z_hat.data())),
            reinterpret_cast<cufftDoubleComplex*>(thrust::raw_pointer_cast(z.data())),
            CUFFT_INVERSE);
        #else        
        cufftExecC2C(plan,
            reinterpret_cast<cufftComplex*>(thrust::raw_pointer_cast(z_hat.data())),
            reinterpret_cast<cufftComplex*>(thrust::raw_pointer_cast(z.data())),
            CUFFT_INVERSE);
        #endif

        normalize<<<(h_N+255)/256, 256>>>(thrust::raw_pointer_cast(z.data()));
        
        REAL Dt=dt;
        thrust::copy(z.begin(), z.end(), z_prev.begin());
        thrust::transform(z.begin(), z.end(), z_prev.begin(), dzdt.begin(),
            [Dt]__device__ __host__ (complex z2, complex z1) {
                return (z2 - z1)/Dt; // dz/dt = (z - z_prev) / dt
            });
    }

    thrust::tuple<complex,complex> rough()
    {
      return roughness(z);
    }
    
    void update_histogram()
    {
      complex one_part_sol = one_particle_solution(h_Ba);
      REAL u0 = one_part_sol.real();
      REAL phi0 = one_part_sol.imag();
      //histogramKernel(z, int* bins, int N, int Nbins, float xmin, float xmax);
    }

    void print_Sq_vs_t(std::ofstream &out, const REAL t)
    {
      complex unit = complex(0.0f, 1.0f);
      
      for(int i = 0; i < h_N/2; ++i) {
        complex z_hat_i = z_hat[i];
        complex z_hat_i_neg = z_hat[h_N-i];
        
        // Note: z_hat[N-i] is the negative frequency component
        complex z_hat_i_neg_conj = thrust::conj(z_hat_i_neg);
        
        complex zu = (z_hat_i + z_hat_i_neg_conj) / 2.0f;  
        complex zphi = unit*(z_hat_i - z_hat_i_neg_conj) / 2.0f;
        
        REAL Squ = zu.real() * zu.real() + zu.imag() * zu.imag();
        REAL Sqphi = zphi.real() * zphi.real() + zphi.imag() * zphi.imag();
        out << 2*M_PI*i/L << " " << Squ << " " << Sqphi << " " << t << std::endl;
      }
      out << "\n" << std::endl;
    }

    cufftHandle plan;
    thrust::device_vector<complex> z;
    thrust::device_vector<complex> z_prev;
    thrust::device_vector<complex> dzdt;
    thrust::device_vector<complex> z_hat;
    thrust::device_vector<complex> nonlinear;
    thrust::device_vector<complex> L_k;
    thrust::device_vector<complex> zaux;
    thrust::device_vector<int> histogram_udot;
    thrust::device_vector<int> histogram_phidot;    
};

int two_system()
{
    Cuerda cuerda1;
    Cuerda cuerda2;

    std::ofstream outz("averagedistances.dat");

    int measurements=0;
    int stride = 1000; // Number of steps between measurements
    thrust::host_vector<complex> distances(stride,complex(0.0f,0.0f));    

    // equilibration
    int eq_steps = 50000;
    for (int n = 0; n < eq_steps; ++n) cuerda1.step();

    // Lyapunov measurements (stride must divide steps)
    for (int n = 0; n < steps; ++n) {
        if(n%stride==0){
          cuerda2.copy_conf(cuerda1);
          cuerda2.perturb_conf(0.0001);
          measurements++;
        }

        cuerda1.step();
        cuerda2.step();

        REAL dist = cuerda1.distance_conf(cuerda2);
        REAL dist_fourier = cuerda1.distance_conf_fourier(cuerda2,16);
        std::cout << dist << " " << dist_fourier << std::endl;
        
        distances[n % stride] += complex(dist, dist_fourier); 
    }
    
    for (int i = 0; i < stride; ++i) {
        outz 
        << distances[i].real()/measurements << " " 
        << distances[i].imag()/measurements << 
        std::endl;
    }
    
    return 0;
}


int one_system()
{
    Cuerda cuerda;
    //cuerda.init();

    std::ofstream outz("z_vs_t.dat");
    std::ofstream outSq("Sq_vs_t.dat");

    long Nmes=0;
    complex av_cm=complex(0.0f, 0.0f);
    complex av_cm2=complex(0.0f, 0.0f);
    complex zcm_middle(0.0f, 0.0f);
    complex zcm_half_steps(0.0f, 0.0f);
    long unsigned int n_middle = 0, n_half_steps = 0;

    unsigned int nlog = 1;


    complex zcm2, zcm;
    thrust::tuple<complex,complex> result = cuerda.rough();
    zcm2 = thrust::get<0>(result);
    zcm = thrust::get<1>(result);
    bool equilibrated = false;
    
    int n = 0;
    for (n = 0; ((zcm.imag()<2.0f*M_PI*100.0f) && n<steps); ++n) {
        cuerda.step();

        if (n % 100 == 0) {
            result = cuerda.rough();
            zcm2 = thrust::get<0>(result);
            zcm = thrust::get<1>(result);
            outz
            << n*dt << " "
            << zcm.real() << " "
            << zcm.imag() << " "
            << zcm2.real() << " "
            << zcm2.imag() << std::endl;
            
            if(zcm.imag()>2*M_PI*5)
            {
              av_cm+=zcm;
              av_cm2+=zcm2;
              Nmes++;
            }
        }
        
        if(zcm.imag()>2*M_PI*90.0f && !equilibrated)
        {
            std::cout << "Equilibration finished at step: " << n << ", distance traveled: " << zcm.imag() << std::endl;
            equilibrated = true;
            
            result = cuerda.rough();
            zcm_middle = thrust::get<1>(result);
            n_middle = n;
        }
        if(n==int(steps*0.8)) {
            result = cuerda.rough();
            zcm_half_steps = thrust::get<1>(result);
            n_half_steps = n;        
        }
        
        if(n == nlog) 
        {
            //outSq << "Step: " << n <<  std::endl;
            nlog *= 2;
            cuerda.print_Sq_vs_t(outSq, n*dt);
        }
    }

    // Save final result
    thrust::host_vector<complex> z_final = cuerda.z;
    thrust::host_vector<complex> dzdt_final = cuerda.dzdt;
    std::ofstream out("z_final.dat");
    for (int i = 0; i < h_N; ++i) {
        out 
        << i * dx << "\t" << z_final[i].real() << "\t" << z_final[i].imag() 
        << "\t" << dzdt_final[i].real() << "\t" << dzdt_final[i].imag()
        << "\n";
    }
    out.close();

    result = cuerda.rough();
    zcm2 = thrust::get<0>(result);
    zcm = thrust::get<1>(result);
    av_cm*=1.0f/Nmes;
    av_cm2*=1.0f/Nmes;

    std::cout << "Run finished at step: " << n << ", distance traveled: " << zcm.imag() << std::endl;

    std::ofstream out_av("averages.dat");

    complex one_part_sol = one_particle_solution(h_Ba);
    
    complex delta_zcm;
    if(equilibrated)
    delta_zcm = (zcm - zcm_middle)/((n-n_middle)*dt);
    else 
    delta_zcm = (zcm-zcm_half_steps)/((n-n_half_steps)*dt);

    out_av
    << h_Ba
    << " " << -delta_zcm.real() << " " << delta_zcm.imag()
    << " " << av_cm2.real() << " " << av_cm2.imag()
    << " " << zcm2.real() << " " << zcm2.imag()
    << " " << one_part_sol.real() << " " << one_part_sol.imag()
    << std::endl;

    /*std::cout
    << h_Ba
    << " " << -delta_zcm.real() << " " << delta_zcm.imag()
    << " " << av_cm2.real() << " " << av_cm2.imag()
    << " " << zcm2.real() << " " << zcm2.imag()
    << " " << one_part_sol.real() << " " << one_part_sol.imag()
    << std::endl;
    */
    
    return 0;
}


int main(int argc, char **argv) {

    // Copy to device constant memory
    h_N = atoi(argv[1]);
    cudaMemcpyToSymbol(N, &h_N, sizeof(int));
    
    // Copy to device constant memory
    h_Ba = atof(argv[2]);
    cudaMemcpyToSymbol(B_a, &h_Ba, sizeof(REAL));
    
    if(argc > 3)
    steps = atoi(argv[3]);
                
    L = h_N*1.0f;
    dx = L / h_N;  
    dt = 3.0f;
    //steps = 500000;

    alpha=complex(0.27f, 0.0f);
    K = 0.796f;
    N_n = 0.016f;

    dt = (h_Ba>alpha.real()*N_n/2.0)?
         (0.01/h_Ba):(0.01/(alpha.real()*N_n/2.0));

    assert(h_Ba*dt < 0.05);
    assert(alpha.real()*N_n*dt/2.0 < 0.05);

    std::ofstream out("parameters.dat");
    out << "N=" << h_N << std::endl;
    out << "B_a=" << h_Ba << std::endl;   
    out << "L=" << L << std::endl;
    out << "dx=" << dx << std::endl;
    out << "dt=" << dt << std::endl;
    out << "steps=" << steps << std::endl;
    out << "alpha=" << alpha << std::endl;
    out << "K=" << K << std::endl;
    out << "N_n=" << N_n << std::endl;
    out << "Using " << (sizeof(REAL) == sizeof(double) ? "double" : "float") << " precision." << std::endl;
    out << "EPSILON=" << EPSILON << std::endl;
    out.close();
    
    #ifndef TWO_SYSTEMS
    one_system();
    #else
    two_system();
    #endif

    return 0;
}
