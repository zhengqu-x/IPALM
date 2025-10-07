# IPALM - Inexact Proximal Augmented Lagrangian Framework

This repository contains the implementation and experimental code for the paper:

**"An inexact proximal augmented Lagrangian framework with arbitrary linearly convergent inner solver for composite convex optimization"** by F. Li and Z. Qu.

## Overview

The IPALM package provides algorithms to solve composite convex optimization problems of the form:

```
minimize f(x) + g(x) + h(p(x))
```

where:
- `f(x)` is convex and differentiable on an open set containing dom(g)
- `g(x)` is proper, convex and closed
- `h(x)` is proper, convex and closed or an indicator function of a convex and closed set
- `p(x)` is a linear operator (matrix) or quadratic operator

## Problem Types Supported

1. **Basis pursuit**: `f(x) = 0`, `g(x) = λ₁‖x‖₁`, `h(x) = I_{x=b}`
2. **Least absolute deviation**: `f(x) = 0`, `g(x) = λ₁‖x‖₁`, `h(x) = λ₃‖x‖₁`
3. **Fused lasso**: `f(x) = 0.5‖Ax - b‖₂²`, `g(x) = λ₁‖x‖₁`, `h(x) = λ₃‖x‖₁`
4. **L₁ norm support vector machine**: `f(x) = 0`, `g(x,w) = λ₁‖x‖₁`, `h(x) = max(0,1-x)`
5. **Quadratically constrained quadratic program**: `f(x) = 0.5xᵀQ₀x + b₀x`, `g(x) = I_{-b ≤ x ≤ b}`, `h(x) = I_{x ≤ 0}`, `p(x) = (0.5xᵀQⱼx + bⱼx - 1)ⱼ₌₁,...,ₘ`

## Repository Structure

### 📁 Core Implementation
- **`IPALM/`** - C++ implementation of five optimization algorithms:
  - IPALM_APPROX
  - IPALM_Katyusha  
  - SMART_CD
  - ASGARD_DL
  - LADMM

### 📁 Data and Benchmarks
- **`datas/`** - 20 benchmark datasets including:
  - Standard datasets: news20scale2, rcv1, rcv2, rcv1mc, news20binary, ijcnn1, w4a-w8a, a7a-a9a, covtype
  - Generated datasets: qcqp1-qcqp4 (quadratically constrained quadratic programs)
  - Some datasets require unzipping or separate download (real-sim, covtype)

### 📁 Analysis and Visualization
- **`MATLAB_code/`** - MATLAB scripts to generate performance plots from experimental results
- **`CVX/`** - CVX package implementations for comparison
- **`libsvm/`** - LIBSVM package implementations for comparison

### 📋 Documentation
- **`EXPERIMENTS.md`** - Detailed instructions for reproducing all paper results and figures

## Quick Start

### Prerequisites
- C++ compiler with GSL (GNU Scientific Library) support
- MATLAB (for generating plots)
- CVX package (optional, for comparison)
- LIBSVM package (optional, for comparison)

### Installation and Compilation

```bash
cd IPALM
g++ -o main main.cpp -lgsl -lgslcblas
```

### Basic Usage

```bash
./main arg1 arg2 arg3 arg4 arg5 arg6 arg7 arg8 arg9 arg10
```

**Parameters:**
- `arg1`: Problem type (1-5)
- `arg2`: Algorithm choice (a-e)
- `arg3`: Dataset filename
- `arg4`: λ₁ parameter
- `arg5`: λ₃ parameter  
- `arg6`: Maximum outer iterations
- `arg7`: Maximum running time (seconds)
- `arg8`: β₀ parameter (optional, default=1)
- `arg9`: ε₀ parameter (optional, default=0)
- `arg10`: Print interval (optional, default=100)

### Example

Solve least absolute deviation problem with IPALM_APPROX:

```bash
./main 2 a news20scale2 0.01 1 100000 100 1 100 100
```

Solve the same problem with SMART_CD:

```bash
./main 2 c news20scale2 0.01 1 100000 100 1000 50
```

## Algorithm Overview

| Code | Algorithm | Description |
|------|-----------|-------------|
| `a` | IPALM_APPROX | Inexact proximal augmented Lagrangian with approximation |
| `b` | IPALM_Katyusha | IPALM with Katyusha inner solver |
| `c` | SMART_CD | Accelerated coordinate descent method |
| `d` | ASGARD_DL | Adaptive primal-dual framework |
| `e` | LADMM | Linearized alternating direction method of multipliers |

## Reproducing Experimental Results

The repository includes detailed instructions for reproducing all figures and tables from the paper. Each experiment involves:

1. Running the C++ executable with specific parameters
2. Processing results with MATLAB visualization scripts

📋 **For complete experimental details and step-by-step reproduction instructions, see [EXPERIMENTS.md](EXPERIMENTS.md)**

### Example: Reproducing Figure 1(a)

```bash
cd IPALM
g++ -o main main.cpp -lgsl -lgslcblas

# Run experiments (takes ~4400 seconds total)
./main 2 a news20scale2 0.01 1 100000 800 1 100 100
./main 2 c news20scale2 0.01 1 100000 1200 1000 50
./main 2 d news20scale2 0.01 1 100000 1200 10
./main 2 e news20scale2 0.01 1 1000000 1200 0.1 100

# Generate plot
cd ../MATLAB_code
matlab -r "plot_news20scale2_lasso"
```

### Results Storage
- Experimental results are saved in `IPALM/results/`
- Generated plots are saved as EPS files in `MATLAB_code/myplots/`

## Dataset Information

### Standard Datasets
Most datasets are sourced from the [LIBSVM repository](https://www.csie.ntu.edu.tw/~cjlin/libsvm/):
- Classification: a7a-a9a, w4a-w8a, ijcnn1, covtype
- Text processing: news20scale2, news20binary, rcv1, rcv1mc

### Generated Datasets  
- `qcqp1`, `qcqp2`: Included in repository
- `qcqp3`, `qcqp4`: Generate using `CVX/qcqp_generate.m` (too large for upload)
- `real-sim`: Download separately from LIBSVM (too large for upload)

### Preprocessing Requirements
Some datasets require extraction:
```bash
# Unzip compressed datasets
unzip matrix_news20binary.zip
unzip matrix_news20binary_bp.zip  
unzip matrix_covtype.zip
```

## Performance Benchmarks

The package includes comprehensive performance comparisons across:
- 5 optimization algorithms
- 20+ benchmark datasets
- Multiple problem formulations (basis pursuit, fused lasso, SVM, QCQP)

Typical runtime expectations:
- Small problems (w4a, w6a): ~200-400 seconds per algorithm
- Medium problems (news20scale2, rcv1): ~600-2400 seconds per algorithm  
- Large problems (news20binary, real-sim): ~3600-14800 seconds per algorithm

## Citation

If you use this code in your research, please cite:

```bibtex
@article{li2021inexact,
  title={An inexact proximal augmented Lagrangian framework with arbitrary linearly convergent inner solver for composite convex optimization},
  author={Li, F. and Qu, Z.},
  journal={Mathematical Programming Computation},
  volume={13},
  pages={583--644},
  year={2021},
  doi={10.1007/s12532-021-00205-x}
}
```

## Dependencies

- **Required**: GSL (GNU Scientific Library), C++ compiler
- **Optional**: MATLAB, CVX package, LIBSVM

## License

This software is licensed for **academic research use only**. See below for details:

```
Academic Research License

Copyright (c) 2021 F. Li and Z. Qu

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to use
the Software solely for academic research, educational, and non-commercial
purposes, subject to the following conditions:

1. The above copyright notice and this permission notice shall be included in 
   all copies or substantial portions of the Software.

2. The Software may only be used for academic research, educational purposes,
   and other non-commercial activities. Commercial use, including but not 
   limited to integration into commercial products or services, is strictly
   prohibited without explicit written permission from the authors.

3. Any publications, presentations, or other academic works that result from
   the use of this Software must cite the original paper:
   Li, F., Qu, Z. An inexact proximal augmented Lagrangian framework with 
   arbitrary linearly convergent inner solver for composite convex optimization.
   Math. Prog. Comp. 13, 583–644 (2021).

4. Redistributions of the Software must retain the above copyright notice,
   this list of conditions, and the following disclaimer.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

For commercial use or licensing inquiries, please contact the authors.
```

## Contact

For questions or issues, please contact the authors or open an issue in this repository.
