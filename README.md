Structure of the Project
=
This project provides subroutines to compute the Jordan-Chevalley decomposition of integer matrices. The implementation consists of a C++ backend a Python frontend.

C++ Backend 
- 
The backend relies heavily on the C++ libraries [Eigen](http://eigen.tuxfamily.org/) and [Boost.Multiprecision](https://www.boost.org/doc/libs/1_72_0/libs/multiprecision/doc/html/index.html). Eigen is a popular algebra library. Boost.Multiprecision is used for variable precision computations and interacts nicely with Eigen. Both are templated header-only libraries and most of their complexity is abstracted away at compile time. 

This project contains a small header-only library which provides simple subroutines as well as a selection of more advanced algorithms for polynomial manipulations. It is used to implement the computing the Jordan-Chevalley decomposition.

Python Frontend
-
The fronted is built using [Cython](https://cython.org) and [Eigency](https://github.com/wouterboomsma/eigency). Cython simplifies the process of writing interfaces between Python and C++ by taking care of the low-level implementation. Eigency uses Cython to implement an interface between between NumPy and Eigen datatypes. The interface is first compiled to C++ and into a shared object file using the Python/C API. 


Setup Instructions
=

I tested the setup on 

* Ubuntu 18.04 with GCC 7.5 (and clang 6.0 for the backend)
* CentOS 7.5 with GCC 6.3
* Windows 10 version 2004 with MSVC 19.25
* macOS 10.15 with Apple clang 11.0.3

Start by changing the current working directory to the directory containing this document.

Setting up Conda Environment
-

An easy way to manage Python environments is provided by Conda. A minimal installation of Conda is provided by [Miniconda](https://docs.conda.io/en/latest/miniconda.html). Please
follow the installation instructions. All default settings are reasonable, except auto_activate_base, deactivate it with:
```
conda config --set auto_activate_base false
```

Now, set up the conda environment


1. Update conda and all packages
    ```
    conda update conda
    conda update --all
    ```

2. Create a conda environment with the packages specified in conda-env.yml.

    I recommend to store the environment in the project folder
      ```
          conda env create -p ./conda-env --file conda-env.yml
      ```  
      and to set the promt to (conda):
      ```
          conda config --set env_prompt '(conda)'
      ``` 

3. Activate the virtual environment with
    ```
    conda activate ./conda-env 
    ```


4. Since eigency is not in the conda repositories, we need to install it manually 
    ```
    pip install eigency
    ```

5. If you wish to deactivate the enviroment, use
    ```
    conda deactivate
    ```

Find more information on managing conda environments [here](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html).



Installing Eigen and Boost
-
There are many ways to install Eigen and Boost, but Conda provides an easy way to install the latest versions of both libraries. Activate the Conda environment, then run

```
conda install -c conda-forge eigen boost-cpp
```
Note that Boost >=1.68 required.


Installing the Python Frontend
-
Assuming that a C++ compiler is available, activate the Conda environment and install with
```
pip install libjordanchevalley/
```
This will compile and install the python extension in the installation into the environment directory and takes roughly 2 minutes.

Uninstalling
-

If you want to just the extension, use
```
pip uninstall jordanchevalley
```


Everything is stored in the Conda environment. Remove it with
```
conda remove -p ./conda-env --all
```

Then, follow the [instructions](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html#uninstalling-anaconda-or-miniconda) to uninstall Miniconda.
 


Using the Extension
=

Please refer to `demo.py` for examples on how to use the extension.