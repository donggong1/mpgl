# mpgl
An *efficient* implementation for the MPGL method (a matching pursuit method for [generalized LASSO](https://arxiv.org/pdf/1005.1971) problem) in the paper:  
<small>*MPGL: An Efficient Matching Pursuit Method for Generalized LASSO  
Dong Gong, Mingkui Tan, Yanning Zhang, Anton van den Hengel, Qinfeng Shi.  
In Thirty-First AAAI Conference on Artificial Intelligence (AAAI), 2017.*  
\[[Paper](https://donggong1.github.io/docs/mpgl_aaai17.pdf)\]\[[Project](https://donggong1.github.io/mpgl.html)\]
</small>

+ If you use this code for your  research, please cite our paper:
````
@inproceedings{gong2017mpgl,
  title={MPGL: An Efficient Matching Pursuit Method for Generalized LASSO},
  author={Gong, Dong and Tan, Mingkui and Zhang, Yanning and van den Hengel, Anton and Shi, Qinfeng},
  booktitle={AAAI Conference on Artificial Intelligence},
  year={2017}
}
````

+ This implementation is based on MATLAB and C++ with an MEX interface. The main framework and the entry of the algorithm is implemented in MATLAB. The most time consuming core part is implemented in C++ with MEX interface for high efficiency.

+ As shown in our paper, the MPGL algorithm and this implementation is both efficient and effective for many applications related to the generalized LASSO problem. 


## Usage
### Getting Started
#### Install [OpenBLAS](http://www.openblas.net/)
+ Download the [source](https://github.com/xianyi/OpenBLAS) to *your/openblas/path/*;
+ Run `Make`;
+ Copy the required *.h* files (i.e. *cblas.h* and *openblas_config.h*) to the *mpgl/mxsolver/*.
+ More details for compiling and installation of OpenBLAS can be found from https://github.com/xianyi/OpenBLAS.


#### Compile *MPGL* solver in C++
+ Set `openblas_path = 'your/openblas/path/'` in *makemexfiles.m*;
+ Run `makemexfiles` in MATLAB and get the *.mex** files.


#### Test MPGL
+ Run `script_gendata` to generate synthetic data for testing. 
+ The entry of the MPGL method is in *entry_mpgl.m*.
+ Please find the instructions for tuning the parameters for data generation and MPGL in the comments in the *.m* files.

#### Solvers
+ This implementation mainly focuses on the model `y=Ax+n` and the generalized LASSO problem. 
+ The solver for the general case of A is in the *mpgl/solver_ls/*.
+ For the cases with `A=I` (i.e. A is an identity matrix, FLSA problem or projection problem), we specifically provide an implementation to accelerate the computation. The related solver is in *mpgl/solver_proj/*.
+ The core parts of the solver in C++ are in *mpgl/mexsolver/*.


### Related Projects
To be updated...