# Installation

### Install Bazel
This gensas repository uses [Bazel](https://bazel.build/) as the build tool. The recommended way to install Bazel is [Bazelisk](https://github.com/bazelbuild/bazelisk). Follow the steps to install bazelisk.

There are also [other ways to install Bazel](https://bazel.build/start).

[//]: # (TODO - rygx: build on other OS (Windows and MacOS) to be verified as a part of https://github.com/ANL-CEEESA/gensas/issues/2)
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
[//]: # (TODO - https://github.com/ANL-CEEESA/gensas/issues/4: fix HDF5 runtime issue)
This is still working in progress (WIP). 