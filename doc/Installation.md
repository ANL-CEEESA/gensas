# Installation

### Install Bazel
This gensas repository uses [Bazel](https://bazel.build/) as the build tool. The recommended way to install Bazel is [Bazelisk](https://github.com/bazelbuild/bazelisk). Follow the steps to install bazelisk.

There are also [other ways to install Bazel](https://bazel.build/start).

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