# Installation

### Dependencies
You will need to install the following dependencies first:

1. [SuperLU](https://portal.nersc.gov/project/sparse/superlu/).

2. [Armadillo](http://arma.sourceforge.net/). Armadillo .hpp files need to be included and build together with GenSAS.

3. BLAS and Lapack. For the best performance, [OpenBLAS](https://www.openblas.net/) is recommended.

4. [MatIO](https://github.com/tbeu/matio). 

5. [HDF5](https://www.hdfgroup.org/solutions/hdf5/).

6. [zlib](https://www.zlib.net/).

7. [jsoncpp](https://github.com/open-source-parsers/jsoncpp).

8. [nvwa] (https://github.com/adah1972/nvwa/tree/master/nvwa). Only `_nvwa.h`, `c++11.h`, `debug_new.h`, `fast_mutex.h`, `pctimer.h`, `static_assert.h` are needed to be included.

Then compile and build GenSAS by including or linking the above dependencies.