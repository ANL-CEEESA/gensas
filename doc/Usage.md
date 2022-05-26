# Usage

### PowerSAS
Current version only provides AC power flow. Assume the executable is gensas.out, then the power flow computation can be called as follows.

```bash
gensas.out -p 
    -f/--file <input-file-name> 
    -o/--output <output-file-name> 
    [-l/--level <sas-order>]
    [-s/--segment <segment-length>]
    [-a/--alphatol <alpha-tolerance>]
    [-d/--difftol <error-tolerance>]
```

Explanations:
* The first argument `-p` (mandatory, and must be the first argument) means GenSAS runs under PowerSAS mode. For the rest of the arguments, the order does not matter.
* `-f/--file <input-file-name>` (mandatory) specifies the input data. Currently the input file needs to be .mat file containing PSAT data structure. See `resources/psat_mat/d_014_syn_ind_zip_export.mat` for example.
* `-o/--output <output-file-name>` (mandatory) specifies the output curve data. Currently the output is a .mat file containing the power flow solution as a vector.
* `-l/--level <sas-order>` (optional) specifies the order of SAS. If not specified, the order of SAS is 15.
* `-s/--segment <segment-length>` (optional) specfies the length of a segment of SAS computation. If not specified, the the segment length is 1.0.
* `-a/--alphatol <alpha-tolerance>` (optional) specifies the tolerance of embedding variable in SAS. If not specified, the tolerance is set as 1e-4.
* `-d/--difftol <error-tolerance>` (optional) specifies the error tolerance of the equations. If not specified, the error tolerance is set as 1e-6.

Example:
Try running power flow of the modified synthetic eastern-interconnection (EI) 70,000-bus system:
```bash
./gensas.out -p -f resources/psat_mat/d_70k_070.mat -s 0.5 -l 28 -d 1e-5 -o res.mat
```

### ModelicaSAS
Currently, ModelicaSAS supports simulation of a single Modelica .mo model without discrete events. The simulation can be called as follows:

```bash
gensas.out -g 
    -m/--mode file/string 
    -i/--input <input> 
    -o/--output <output-file-name>  
    [-j/--json <json-output-file-name>]
    [-l/--level <sas-order>]
    [-s/--segment <segment-length>]
    [-a/--aTol <alpha-tolerance>]
    [-e/--eTol <error-tolerance>]
    [-t/--outStep <output-step>]
    [-v/--verbose]
```

Explanations:
* The first argument `-g` (mandatory, and must be the first argument) means GenSAS runs under ModelicaSAS mode. For the rest of the arguments, the order does not matter.
* `-m/--mode file/string` (mandatory) specifies the format of input. If `-m file`, then expect the input file name after `-i`; and if `-m string`, then expect the Modelica model content after `-i`.
* `-i/--input <input>` (mandatory) specifies the input data. The `<input>` depends on the mode specified after `-m`.
* `-o/--output <output-file-name>` (mandatory) specifies the output curve data file name. Currently the output is a .mat file containing the output as a matrix.
* `-j/--json <json-output-file-name>` (optional) specifies the name of the .json file containing output curve. 
* `-l/--level <sas-order>` (optional) specifies the order of SAS. If not specified, the order of SAS is 15.
* `-t/--time <max-time>` (optional) specifies the maximum time of simulated process. If not specified, the time is set as 10.0.
* `-s/--segment <segment-length>` (optional) specfies the length of a segment of SAS computation. If not specified, the the segment length is 1.0.
* `-a/--atol <alpha-tolerance>` (optional) specifies the tolerance of embedding variable in SAS. If not specified, the tolerance is set as 1e-3.
* `-e/--etol <error-tolerance>` (optional) specifies the error tolerance of the equations. If not specified, the error tolerance is set as 1e-5.
* `-k/--step <output-step>` (optional) specifies the time step of the output curves. If not specified, the time step is set as 0.01.
* `-v/--verbose` (optional) if used, will print intermediate result in SAS computation.

Example:
Try running simulation of the model in `resources/mofile/test_solve_ode.mo`.

```bash
 ./gensas.out -g 
    -m file -i resources/mofile/test_solve_ode.mo 
    -o resources/mofile/test_solve_ode.mat 
    -t 15
```