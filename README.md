# MDOODZ6.0

Welcome to MDOODZ6.0 public repository!

# Quickstart

To reproduce the results of the (submitted) manuscript:

Check the Documentation

Compile: make clean all MODEL=LithoScale OPT=yes OMP=yes

Run: ./Doodzi_LithoScale LithoScale.txt

Visualize using either Matlab, Python, Julia or whatever language that can handle HDF5 files and enjoy!

# Prequisites

So far MDOODZ has been successfully built on LINUX/UNIX and MAC OS systems. The code can be built with GCC compiler from GNU (http://gcc.gnu.org) or with ICC compiler from Intel MKL library (https://software.intel.com/en-us/intel-mkl).
The code relies on two libraries: <br>
1. SuiteSparse provides efficient linear algebra and matrix manipulation routines. It is available at: http://www.suitesparse.com <br>
2. HDF5 is the main format for output files and is readable into MATLAB. It is available at: http://www.hdfgroup.org <br>


