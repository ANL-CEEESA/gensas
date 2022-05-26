# GenSAS

**GenSAS** is a C++ based generalized simulation tool based on semi-analytical solutions (SAS) technology. It provides two modules:

1. **ModelicaSAS** is a SAS-based simulation tool that supports Modelica modeling language. Currently it supports simulation of single Modelica model.

2. **PowerSAS** is a robust, efficient and scalable power grid analysis framework based on semi-analytical solutions (SAS) technology. This is the version for C++ users. It currently provides AC power flow functionality.

Backed by the SAS approach, the GenSAS tool provides much better convergence than the tools using traditional Newton-type equation solvers. Moreover, due to the analytical nature, GenSAS provides model-adaptive high-accuracy approximation, which brings significantly extended effective range for power grid analysis. PowerSAS has been used to solve large-scale system cases with 200,000+ buses.

### Companion Project
PowerSAS.m -- SAS-based power grid analysis toolbox for Matlab/GNU Octave users. 

### Acknowledgement
This work is supported by the Laboratory Directed Research and Development (LDRD) program of Argonne National Laboratory, provided by the U.S. Department of Energy Office of Science under Contract No. DE-AC02-06CH11357, and the U.S. Department of Energy Office of Electricity, Advanced Grid Modeling program under Grant DE-OE0000875.

### Publications
* Rui Yao, Kai Sun, Feng Qiu, “Vectorized Efficient Computation of Padé Approximation for Semi-Analytical Simulation of Large-Scale Power Systems,” IEEE Transactions on Power Systems, 34 (5), 3957-3959, 2019.
* Rui Yao, Yang Liu, Kai Sun, Feng Qiu, Jianhui Wang,"Efficient and Robust Dynamic Simulation of Power Systems with Holomorphic Embedding", IEEE Transactions on Power Systems, 35 (2), 938 - 949, 2020.
* Rui Yao, Feng Qiu, "Novel AC Distribution Factor for Efficient Outage Analysis", IEEE Transactions on Power Systems, 35 (6), 4960-4963, 2020.
* Xin Xu, Rui Yao, Kai Sun, Feng Qiu, "A Semi-Analytical Solution Approach for Solving Constant-Coefficient First-Order Partial Differential Equations", IEEE Control Systems Letters, 6, 704-709, 2021.
* Rui Yao, Feng Qiu, Kai Sun, “Contingency Analysis Based on Partitioned and Parallel Holomorphic Embedding”, IEEE Transactions on Power Systems, in press.

### License
This software is under 3-clause BSD license. Please refer to LICENSE.md for details.