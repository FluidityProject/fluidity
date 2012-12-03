
from numpy import linalg
from numpy import zeros, matrix

lms=matrix(zeros((12,12)))

f=('local_solver_mat1', 'r')
for i in range(12):
    l=f.readline().split()
    for j in range(12):
        lms[i,j]=l[j]

rhs=matrix(zeros((12)))

frhs=('fort.1', 'r')
for i in range(12):
    l=frhs.readline()
