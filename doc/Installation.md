# Installation

GenSAS builds with [Bazel](https://bazel.build/) on Linux, macOS and Windows.

| | Linux | macOS | Windows |
|---|---|---|---|
| Toolchain | GCC or Clang, `gfortran`, `cmake` | Xcode command line tools, Homebrew `gcc` (provides `gfortran`), `cmake` | MSYS2 MinGW-w64 GCC, CMake, Ninja (no Fortran needed) |
| Shell for the commands below | any | any | **MSYS2 MINGW64** |
| Binary after build | `bazel-bin/app/app` | `bazel-bin/app/app` | `bazel-bin/app/app.exe` |
| Status | built and tested in CI | built and tested in CI | built in CI (non-blocking until it has a green run); verified locally with MSYS2 GCC 16.2 |

### Install Bazel
This gensas repository uses [Bazel](https://bazel.build/) as the build tool. The recommended way to install Bazel is [Bazelisk](https://github.com/bazelbuild/bazelisk). Follow the steps to install bazelisk. The repository pins the Bazel version in `.bazelversion`.

There are also [other ways to install Bazel](https://bazel.build/start).

### Install prerequisites

* **Linux** (Debian/Ubuntu example): `sudo apt install build-essential gfortran cmake`. On Fedora/RHEL: `sudo dnf install gcc-c++ gcc-gfortran cmake make`.
* **macOS** (Intel and Apple Silicon): `xcode-select --install`, then `brew install gcc cmake`. `.bazelrc` adds the Homebrew directories (`/opt/homebrew` and `/usr/local`) to the build `PATH`.
* **Windows**: see [Windows](#windows) below.

### Build Gensas
At the root of the repository, execute command (on Windows use `bazel build //app:app` from the MSYS2 MINGW64 shell, see [Windows](#windows)):
```shell
bazel build //...
```

The first-time build (or fresh re-build after running `bazel clean`) will take a while to finish.

### Running the application
Run the built binary from the repository root so that relative paths work the same way on every OS:

```shell
bazel-bin/app/app <arguments>        # Linux and macOS
bazel-bin/app/app.exe <arguments>    # Windows (MSYS2 shell)
```

Notes:

* `bazel run //app:app -- <arguments>` also works on Linux and macOS, but the program then starts in Bazel's runfiles directory, so input files need absolute paths (for example `$(pwd)/resources/...`) and a relative `-o` file is written relative to that directory instead of your repository root. On macOS the Bazel sandbox may also block writing to arbitrary locations. The examples below therefore run the binary directly.
* For power flow (`-p`), the `-o` path is relative to the current directory; do not pass an absolute path.
* On Windows, native programs cannot open MSYS-style paths such as `/d/work/gensas/...`. Use relative paths from the repository root, or convert absolute ones with `cygpath -m <path>` (for example `-f "$(cygpath -m "$PWD")/resources/..."`).

### Run Modelica example
To run simulation on a small dynamic system under `resources/mofile/test_solve_ode.mo`, run command at repository root (use `app.exe` on Windows):
```shell
bazel-bin/app/app -g -m file -i resources/mofile/test_solve_ode.mo -j solution.json
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

An example for running the synthetic 210,000-bus system (use `app.exe` on Windows):
``` shell
bazel-bin/app/app -p -f resources/psat_mat/d_210k.mat -l 25 -s 0.2 -o res.mat
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
`C_LAPACK=ON`, so **no Fortran compiler is needed on Windows**. OpenBLAS is built
for the `HASWELL` target (AVX2), so a CPU with AVX2 is required.

1. Install [MSYS2](https://www.msys2.org/) to `C:\msys64` and, from an
   "MSYS2 MINGW64" shell, the toolchain:
   ```shell
   pacman -S --needed mingw-w64-x86_64-gcc mingw-w64-x86_64-cmake \
                      mingw-w64-x86_64-ninja mingw-w64-x86_64-make make
   ```
   Only the MINGW64 environment is supported; the paths in `.bazelrc`
   (`build:windows`) point to `C:/msys64/mingw64`. The MSYS2 installer marks
   MINGW64 as deprecated in favour of UCRT64, but UCRT64 is not tested here.
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

### Local overrides
Settings that are specific to your machine (on any OS), such as a different
MSYS2 location, a disk cache, or extra `--action_env=PATH=...` entries, go into
`.bazelrc.local` at the repository root. It is imported by `.bazelrc` and is
not tracked by git.

### Troubleshooting

* **Linux**: `cannot find -lgfortran` or `gfortran: not found` means the Fortran compiler is missing; install `gfortran` (see the prerequisites above).
* **macOS**: if `gfortran` or `cmake` is not found during the build, check that Homebrew is installed in `/opt/homebrew` or `/usr/local`; for another location, set `--action_env=PATH=...` and `--host_action_env=PATH=...` in `.bazelrc.local`.
* **Windows**: CMake errors about `cmake_minimum_required`: the `.bazelrc` sets
  `CMAKE_POLICY_VERSION_MINIMUM=3.5` for this; make sure the Windows section is
  not overridden.
* **Windows**: `app.exe` fails to start with missing DLL errors: add
  `C:\msys64\mingw64\bin` to `PATH`, or build with `--linkopt=-static`.
* **Windows**: input file not found although the file exists: pass a Windows-style or relative path (see "Running the application").
* **All**: if a build fails while downloading a dependency (HTTP 429/5xx), re-run the build; downloads are cached by Bazel.
