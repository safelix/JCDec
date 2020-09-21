# cython: language_level=3
from eigency.core cimport *
from libcpp cimport bool
from libcpp.string cimport string
cimport numpy as np

# Declaring a C++ class interface
# https://cython.readthedocs.io/en/latest/src/userguide/wrapping_CPlusPlus.html#declaring-a-c-class-interface
cdef extern from "jordanchevalley.hpp":

    cdef cppclass _JCDec "JCDec":
        int verbosity

        # Constructor
        _JCDec( Map[MatrixXd] &, bool, int) except+

        # Attribute accessors
        int get_size() except+
        int get_no_iter() except+
        MatrixXd get_A() except+
        VectorXd get_chiA() except+
        VectorXd get_muD() except+
        VectorXd get_inv() except+
        VectorXd get_chev() except+
        MatrixXd get_D() except+
        MatrixXd get_N() except+

        # Compute methods
        _JCDec& compute_chiA(string, bool) except+
        _JCDec& compute_muD(bool) except+
        _JCDec& compute_inv() except+
        _JCDec& compute_chev() except+
        _JCDec& compute_D(bool) except+
        _JCDec& compute_N() except+
        _JCDec& compute() except+

        # Check methods
        bool is_trivial() except+
        double check_cayleyhamilton() except+
        double check_nillpotency() except+
        double check_commutativity() except+

        # Meta attribute accessors
        double get_timing(string) except+
        RowVectorXd get_timings(int) except+
        double get_maxBitsize(string) except+
        RowVectorXd get_maxBitsizes(int) except+


# Create Cython wrapper class
# https://cython.readthedocs.io/en/latest/src/userguide/wrapping_CPlusPlus.html#create-cython-wrapper-class
cdef class JCDec:
    cdef _JCDec *thisptr;

    # Constructor and Deconstructor
    def __cinit__(self, np.ndarray[np.float64_t, ndim=2] mat, *, int verbosity = 0):
        self.thisptr = new _JCDec( Map[MatrixXd](mat), mat.flags.c_contiguous, verbosity )
    
    def __dealloc__(self):
        del self.thisptr

    # Attribute access
    @property 
    def __size__(self):
        return self.thisptr.get_size()
    @property 
    def __no_iter__(self):
        return self.thisptr.get_no_iter()
    @property
    def __A__(self):
        return ndarray_copy(self.thisptr.get_A())
    @property
    def __chiA__(self):
        return ndarray_copy(self.thisptr.get_chiA())
    @property
    def __muD__(self):
        return ndarray_copy(self.thisptr.get_muD()) 
    @property   
    def __inv__(self):
        return ndarray_copy(self.thisptr.get_inv())
    @property
    def __chev__(self):
        return ndarray_copy(self.thisptr.get_chev())
    @property
    def __D__(self):
        return ndarray_copy(self.thisptr.get_D())
    @property
    def __N__(self):
        return ndarray_copy(self.thisptr.get_N())
    
    @property 
    def __verbosity__(self):
        return self.thisptr.verbosity
    @__verbosity__.setter
    def __verbosity__(self, verbosity):
        self.thisptr.verbosity = verbosity

    # Compute methods
    def compute_chiA(self, *, str precision = "double", bool round=True):
        self.thisptr.compute_chiA(precision.encode('UTF-8'), round)
        return self
    def compute_muD(self, bool resultant=False):
        self.thisptr.compute_muD(resultant)
        return self
    def compute_inv(self):
        self.thisptr.compute_inv()
        return self
    def compute_chev(self):
        self.thisptr.compute_chev()
        return self
    def compute_D(self, *, bool mat=False):
        self.thisptr.compute_D(mat)
        return self
    def compute_N(self):
        self.thisptr.compute_N()
        return self
    def compute(self):
        self.thisptr.compute()
        return self

    # Check methods
    def is_trivial(self):
        return self.thisptr.is_trivial()
    def check_cayleyhamilton(self):
        return self.thisptr.check_cayleyhamilton()
    def check_nillpotency(self):
        return self.thisptr.check_nillpotency()
    def check_commutativity(self):
        return self.thisptr.check_commutativity()

    # Meta attribute accessors
    def get_timing(self, str step_name=""):
        return self.thisptr.get_timing(step_name.encode('UTF-8'))

    def get_timings(self, int steps=7):
        return ndarray_copy(self.thisptr.get_timings(steps))

    def get_maxBitsize(self, str step_name=""):
        return self.thisptr.get_maxBitsize(step_name.encode('UTF-8'))

    def get_maxBitsizes(self, int steps=7):
        return ndarray_copy(self.thisptr.get_maxBitsizes(steps))
