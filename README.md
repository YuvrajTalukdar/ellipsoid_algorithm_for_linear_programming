# ellipsoid_algorithm_for_linear_programming

This is an implementation of the Kachiyan's ellipsoid algorithm using gsl and blas.
The repository also includes blas wrapper code to simplify the blas function calls.
test_func folder contains some code for testing log, gsl and blas code.


To install blas libraries in Ubuntu: sudo apt-get install libblas-dev libblas64-dev libatlas-base-dev liblapack-dev libopenblas-dev libgsl-dev
Other install commands were not working, this one works.

For gls it looks like the libraries were already installed, do verify this in VMs.
