CXX = nvcc


INCLUDES = -I/opt/nvidia/hpc_sdk/Linux_x86_64/23.7/math_libs/12.2/include 
FLAGS = --expt-extended-lambda -lcufft -std=c++17 -lstdc++fs \
-gencode arch=compute_61,code=sm_61 -gencode arch=compute_80,code=sm_80 -gencode arch=compute_75,code=sm_75

PARAMS = -DTWO_SYSTEMS -DDOUBLE_PRECISION 

LDFLAGS = -L/opt/nvidia/hpc_sdk/Linux_x86_64/23.7/math_libs/12.2/lib64 


spatialuphispectral: main.cu
	$(CXX) $(FLAGS) $(PARAMS) main.cu -o spatialuphispectral $(LDFLAGS) $(INCLUDES) 


update_git:
	git add *.cu Makefile *.h *.sh README.md *.gnu; git commit -m "program update"; git push

clean:
	rm spatialuphispectral
