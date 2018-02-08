import numpy as np
import matplotlib.pyplot as plt

#  convert to a column vector
def MakeCol(y): return y.reshape(-1,1)
#  convert to a row vector
def MakeRow(y): return y.reshape(1,-1)

# load a dataset
from scipy import io
name = 'iris'
U = io.loadmat('nt_toolbox/data/ml-' + name)
A = U['A']

# find an element
def find(x): return np.nonzero(x)[0]

### useful fonction
np.linalg.norm(x)
# A*x=y
np.linalg.solve(A,y)
