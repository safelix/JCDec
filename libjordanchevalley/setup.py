#!/usr/bin/env python
from setuptools import setup, Extension
import os, glob
import eigency


__root_dir__ = os.path.dirname(__file__)


__python_dir__ = os.path.join(__root_dir__, "python")
__python_src__ = glob.glob(os.path.join(__python_dir__, "*.pyx"))

__core_dir__ = os.path.join(__root_dir__, "core")
__core_src__ = [os.path.join(__core_dir__, "jordanchevalley.cpp")]

__conda_dir__ =  os.environ['CONDA_PREFIX']
if os.name == "posix":
    __conda_include_dir__ = os.path.join(__conda_dir__, "include")
if os.name == "nt":
    __conda_include_dir__ = os.path.join(__conda_dir__, "Library", "include")
    

__eigen__include_dir__ = os.path.join(__conda_include_dir__, "eigen3")

if os.name == "posix":
   compile_args = ["-std=c++17", "-Wall"]
if os.name == "nt":
    compile_args = ["/std:c++17"] # /Wall is too aggressive

extensions = [
    Extension(
        language="c++",
        name="jordanchevalley",
        sources=__core_src__ + __python_src__,
        include_dirs=[__conda_include_dir__, __eigen__include_dir__, __core_dir__]
        + eigency.get_includes(include_eigen=False),
        extra_compile_args=compile_args,
    )
]

dist = setup(
    name="jordanchevalley",
    install_requires=["numpy"],
    setup_requires=["setuptools>=18.0", "cython", "eigency"],
    ext_modules=extensions,
)

