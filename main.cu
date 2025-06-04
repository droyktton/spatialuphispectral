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

const int N = 8192;
const REAL L = N*1.0f;
const REAL dx = L / N;
const REAL dt = 0.2f;
const int steps = 1000000;

const complex alpha(0.27f, 0.0f);
const REAL K = 0.796f;
const REAL N_n = 0.016f;
REAL h_Ba;
__constant__ REAL B_a;



__global__ void init_wave_numbers(complex* L_k, int N, REAL K, REAL L) {
    int i = threadIdx.x + blockDim.x * blockIdx.x;
    if (i < N) {
        int k = (i <= N/2) ? i : i - N;
        REAL kx = 2 * M_PI * k / L;
        L_k[i] = complex(-K * kx * kx, 0.0f);
    }
}

__global__ void nonlinear_term_kernel(const complex* z, complex* nonlinear, REAL N_n, int N) {
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i < N) {
        REAL phi = -z[i].imag();
        //nonlinear[i] = -complex(0.0f, 1.0f) * (N_n / 2.0f) * sinf(2 * phi); // different in the notes...
        nonlinear[i] = complex(0.0f, 1.0f) * (N_n / 2.0f) * sinf(2 * phi);
    }
}

__global__ void crank_nicholson_update(complex* z_hat, const complex* N_hat,
                                       const complex* L_k,
                                       complex alpha, REAL dt, REAL B_a, int N) {
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

__global__ void normalize(complex* z, int N) {
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i < N) {
        z[i] /= static_cast<REAL>(N);
    }
}

thrust::tuple<complex,complex> roughness(thrust::device_vector<complex> &z)
{
    size_t N = z.size();
    complex zcm = thrust::reduce(z.begin(), z.end());
    zcm *= 1.0f/N;

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

    zcm2 *= 1.0f/N;

    return make_tuple(zcm2, zcm);
}

class Cuerda
{
    public:
    Cuerda(){
      z.resize(N);
      z_hat.resize(N);
      nonlinear.resize(N);
      L_k.resize(N);
      //zaux.resize(N);

      init();

      #ifdef DEBUG
      std::cout << "N=" << N << std::endl;
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

        return dist/N;
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

      return dist/N;
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
      thrust::host_vector<complex> hz(N);
      complex avhz;
      for (int i = 0; i < N; ++i) {
        hz[i] = epsilon*complex(rand()*1.0f/RAND_MAX, rand()*1.0f/RAND_MAX);
        avhz += hz[i];
      }
      avhz *= 1.0f/N;
      for(int i=0; i<N; i++) z[i] += (hz[i]-avhz);
    }

    ~Cuerda(){
      cufftDestroy(plan);
    };

    void init(){
      #ifdef DOUBLE_PRECISION
      cufftPlan1d(&plan, N, CUFFT_Z2Z, 1);
      #else 
      cufftPlan1d(&plan, N, CUFFT_C2C, 1);
      #endif
      
      srand(42);
      // Initial condition: z = cos(x)
      thrust::host_vector<complex> z0(N);
      for (int i = 0; i < N; ++i) {
          REAL x = i * dx;
          //z0[i] = complex(0.0f, 0.0f);
          z0[i] = complex(rand()*1.0f/RAND_MAX, rand()*1.0f/RAND_MAX);
          //if(i<10) std::cout << z0[i] << std::endl;
      }
      z = z0;

      // Init linear operator
      init_wave_numbers<<<(N+255)/256, 256>>>(thrust::raw_pointer_cast(L_k.data()), N, K, L);

      #ifdef DEBUG
      std::cout << "K=" << K << std::endl;
      std::cout << "L=" << L << std::endl;
      std::cout << "dx=" << dx << std::endl;
      std::cout << "dt=" << dt << std::endl;
      std::cout << "steps=" << steps << std::endl;
      std::cout << "alpha=" << alpha << std::endl;
      std::cout << "N_n=" << N_n << std::endl;
      std::cout << "B_a=" << h_Ba << std::endl;
      #endif
    }

    void step(){
       nonlinear_term_kernel<<<(N+255)/256, 256>>>(
            thrust::raw_pointer_cast(z.data()),
            thrust::raw_pointer_cast(nonlinear.data()), N_n, N);

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
      
        crank_nicholson_update<<<(N+255)/256, 256>>>(
            thrust::raw_pointer_cast(z_hat.data()),
            thrust::raw_pointer_cast(nonlinear.data()),
            thrust::raw_pointer_cast(L_k.data()),
            alpha, dt, h_Ba, N);

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

        normalize<<<(N+255)/256, 256>>>(thrust::raw_pointer_cast(z.data()), N);
    }

    thrust::tuple<complex,complex> rough()
    {
      return roughness(z);
    }

    void print_Sq_vs_t(std::ofstream &out, const REAL t)
    {
      for(int i = 0; i < N/2; ++i) {
        complex z_hat_i = z_hat[i];
        REAL q_i = 2*M_PI*i/L;
        REAL Sq = z_hat_i.real() * z_hat_i.real() + z_hat_i.imag() * z_hat_i.imag();
        out << q_i << " " << Sq << " " << " " << t << " " << z_hat_i.real() << " " << z_hat_i.imag() << std::endl;
      }
      out << "\n" << std::endl;
    }

    cufftHandle plan;
    thrust::device_vector<complex> z;
    thrust::device_vector<complex> z_hat;
    thrust::device_vector<complex> nonlinear;
    thrust::device_vector<complex> L_k;
    thrust::device_vector<complex> zaux;

};

int two_system()
{
    Cuerda cuerda1;
    Cuerda cuerda2;

    std::ofstream outz("averagedistances.txt");

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
          cuerda2.perturb_conf(0.001);
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

complex one_particle_solution(REAL h){
    REAL a = alpha.real();
    REAL hw = a*N_n/2.0f; 
    REAL vphi = (h > hw)?(sqrtf((h/hw)*(h/hw)-1.0f)*hw/(1+a*a)):0.0f;  
    REAL vu = (h > hw)?(h/a-sqrtf((h/hw)*(h/hw)-1.0f)*hw/(a+a*a*a)):(h/a);
    complex z(vu, vphi);
    return z;  
}

int one_system()
{
    Cuerda cuerda;
    //cuerda.init();

    std::ofstream outz("z_vs_t.txt");
    std::ofstream outSq("Sq_vs_t.txt");

    long Nmes=0;
    complex av_cm=complex(0.0f, 0.0f);
    complex av_cm2=complex(0.0f, 0.0f);
    complex zcm_middle(0.0f, 0.0f);
    long unsigned int n_middle = int(steps*0.5);

    unsigned int nlog = 1;

    for (int n = 0; n < steps; ++n) {
        cuerda.step();

        if (n % 100 == 0) {
            thrust::tuple<complex,complex> result = cuerda.rough();
            complex zcm2 = thrust::get<0>(result);
            complex zcm = thrust::get<1>(result);
            outz
            << n << " "
            << zcm.real() << " "
            << zcm.imag() << " "
            << zcm2.real() << " "
            << zcm2.imag() << std::endl;
            
            if(n>n_middle)
            {
              av_cm+=zcm;
              av_cm2+=zcm2;
              Nmes++;        
            }
        }
        
        if(n==n_middle)
        {
            thrust::tuple<complex,complex> result = cuerda.rough();
            zcm_middle = thrust::get<1>(result);
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
    std::ofstream out("z_final.txt");
    for (int i = 0; i < N; ++i) {
        out << i * dx << "\t" << z_final[i].real() << "\t" << z_final[i].imag() << "\n";
    }
    out.close();

    thrust::tuple<complex,complex> result = cuerda.rough();
    complex zcm2 = thrust::get<0>(result);
    complex zcm = thrust::get<1>(result);
    av_cm*=1.0f/Nmes;
    av_cm2*=1.0f/Nmes;

    std::ofstream out_av("averages.dat");

    complex one_part_sol = one_particle_solution(h_Ba);
    complex delta_zcm = (zcm - zcm_middle)/((steps-n_middle)*dt);

    out_av
    << h_Ba
    << " " << -delta_zcm.real() << " " << delta_zcm.imag()
    << " " << av_cm2.real() << " " << av_cm2.imag()
    << " " << zcm2.real() << " " << zcm2.imag()
    << " " << one_part_sol.real() << " " << one_part_sol.imag()
    << std::endl;

    std::cout
    << h_Ba
    << " " << -delta_zcm.real() << " " << delta_zcm.imag()
    << " " << av_cm2.real() << " " << av_cm2.imag()
    << " " << zcm2.real() << " " << zcm2.imag()
    << " " << one_part_sol.real() << " " << one_part_sol.imag()
    << std::endl;

    return 0;
}


int main(int agrc, char **argv) {

    // Copy to device constant memory
    h_Ba = atof(argv[1]);
    cudaMemcpyToSymbol(B_a, &h_Ba, sizeof(REAL));

    #ifndef TWO_SYSTEMS
    one_system();
    #else
    two_system();
    #endif

    return 0;
}
