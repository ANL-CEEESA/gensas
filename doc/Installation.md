# Installation

### Install Bazel
This gensas repository uses [Bazel](https://bazel.build/) as the build tool. The recommended way to install Bazel is [Bazelisk](https://github.com/bazelbuild/bazelisk). Follow the steps to install bazelisk.

There are also [other ways to install Bazel](https://bazel.build/start).

### Prerequisites by platform

* **Linux**: a C/C++ toolchain, `gfortran` and `cmake`.
* **macOS**: Xcode command line tools and `brew install gcc cmake`.
* **Windows**: MSYS2 with the MinGW-w64 toolchain, see [Windows](#windows) below.

### Build Gensas
At the root of the repository, execute command:
```shell
bazel build //...
```

The first-time build (or fresh re-build after running `bazel clean`) will take a while to finish.

### Run Modelica example
To run simulation on a small dynamic system under `resources/mofile/test_solve_ode.mo`, run command at repository root:
```shell
bazel run //app:app -- -g -m file -i $(pwd)/resources/mofile/test_solve_ode.mo -j /tmp/solution.json
```

You should see output like
```
compMode=-g
All 2391 characters read successfully.
Number of DE: 3
Number of AE: 0
Number of X: 3
Number of Y: 0
Init values generated.
Stage start 0 end 1, length 1
Stage start 1 end 2, length 1
Stage start 2 end 3, length 1
Stage start 3 end 4, length 1
Stage start 4 end 5, length 1
Stage start 5 end 6, length 1
Stage start 6 end 7, length 1
Stage start 7 end 8, length 1
Stage start 8 end 9, length 1
Stage start 9 end 10, length 1
Computation time: 0.00302792 s.
Json file /tmp/solution.json written.
```

### Run a power flow

Power flow can be run using the `-p` mode. There are sample input data representing some test models, with sizes ranging from 14 buses to 210,000 buses.

An example for running the synthetic 210,000-bus system:
``` shell
bazel run //app:app -- -p -f $(pwd)/resources/psat_mat/d_210k.mat -l 25 -s 0.2
```

You should see output like
```
compMode=-p
Step=0.19375, added=0.19375, (maxDiff<1e-06).
Step=0.39375, added=0.2, (maxDiff<1e-06).
Step=0.553125, added=0.159375, (maxDiff<1e-06).
Step=0.740625, added=0.1875, (maxDiff<1e-06).
Step=0.9, added=0.159375, (maxDiff<1e-06).
Step=1, added=0.1, (maxDiff<1e-06).
Computation No.1
Computation time: 89.9048 s.
```

### Windows

Windows is built natively with the MSYS2 MinGW-w64 GCC toolchain (not MSVC); it
has been verified with MSYS2 GCC 16.2 and Bazel 7.0.0.
Dependencies that use autotools on Linux/macOS (HDF5, libmatio) and `make`
(OpenBLAS) are built with CMake on Windows, and OpenBLAS is built with
`C_LAPACK=ON`, so **no Fortran compiler is needed on Windows**.

1. Install [MSYS2](https://www.msys2.org/) to `C:\msys64` and, from an
   "MSYS2 MINGW64" shell, the toolchain:
   ```shell
   pacman -S --needed mingw-w64-x86_64-gcc mingw-w64-x86_64-cmake \
                      mingw-w64-x86_64-ninja mingw-w64-x86_64-make make
   ```
2. Install Bazel (Bazelisk) for Windows. The repository pins Bazel 7.0.0
   (`.bazelversion`).
3. From the **MSYS2 MINGW64 shell**, at the repository root:
   ```shell
   export MSYS2_ARG_CONV_EXCL="*"   # stop MSYS from rewriting Bazel's arguments
   bazel build //app:app
   ```
   (`rules_foreign_cc` assumes an MSYS2 shell on Windows.) The Windows settings
   in `.bazelrc` (`build:windows`) are applied automatically. If MSYS2 is not in
   `C:\msys64`, copy those lines into `.bazelrc.local` and adjust the paths.
4. Run, for example:
   ```shell
   bazel-bin/app/app.exe -p -f resources/psat_mat/d_014_ind_zip1.mat -s 0.5 -l 28 -d 1e-5 -o out.mat
   ```
   or run `bash scripts/run_integration_test.sh`.

Use `bazel build //app:app` rather than `//...` on Windows: targets that only make
sense with the autotools build (e.g. `//:hdf5_lib`, `//:libaec`, `//:pcre`) are
marked incompatible and are skipped.

Troubleshooting:

* CMake errors about `cmake_minimum_required`: the `.bazelrc` sets
  `CMAKE_POLICY_VERSION_MINIMUM=3.5` for this; make sure the Windows section is
  not overridden.
* `app.exe` fails to start with missing DLL errors: add
  `C:\msys64\mingw64\bin` to `PATH`, or build with `--linkopt=-static`.
