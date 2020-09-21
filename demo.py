import numpy as np
from jordanchevalley import JCDec

# Function Definitions
################################################################################

# generate a Kroneker graph (https://cs.stanford.edu/people/jure/pubs/kronecker-jmlr10.pdf)
def kroneker_graph(k, init):
    P = np.eye(1)
    for _ in range(0,k):
        P = np.kron(init,P)
    return np.random.binomial(n=1, p=P)


# START SCRIPT
################################################################################

# initial kroneker probabilities
init = np.array([[0.9, 0.2],[0.2, 0.9]]) 

# generate Kroneker graph
A = kroneker_graph(2, init)
print(A)

# initiallize Jordan-Chevalley decomposition object
# - set the verbosity to 1
dec = JCDec(A.astype('float64'), verbosity=1)


# compute the characteristic polynomial of A using 
# - 256-bit FP precision (oct)
# - rounding coefficients to next integer
dec.compute_chiA(precision="mp::oct", round=True)


# compute the minimal polynomial of D using
# - the resultant algorithm
dec.compute_muD(resultant=True)


# compute the inv(x) and chev(x)
# - it is possible to create chains of commands
dec.compute_inv().compute_chev()


# compute the diagonalizable matrix D
# - using the Chevalley iteration on matrices
dec.compute_D(mat=True)


# compute the nilpotent matrix N
# - required intermediate results are computed automatically
# dec.compute_N()


# genaral function to compute the decomposition in  one go
dec.compute()


# some checks
if dec.is_trivial():
    print("the input matrix was diagonalizable")
else:
    print("the input matrix was not diagonalizable")


if dec.check_commutativity() == 0 and dec.check_nillpotency() == 0:
    print("the computed decomposition is correct")
else:
    print("the computed decomposition is incorrect")


# inspect (intermediate) results as a numpy matrix
print(dec.__D__)


