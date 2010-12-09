#!/usr/bin/env python

from numpy import *
from numpy.linalg import solve

n = 5

D = zeros((n, n))
D[0,0] = 1
D[0,1] = -1
D[n-1,n-1] = 1
D[n-1,n-2] = -1
for i in range(1,n-1):
  D[i,i] = 2
  D[i,i-1] = -1
  D[i,i+1] = -1

D = (n-1) * D

A = zeros((n, n))
A[0,0] = -1
A[0,1] = 1

A = (n-1) * A

B = zeros((n, n))
B[n-1,n-1] = 1
B[n-1,n-2] = -1

B = (n-1) * B

j = zeros(n)
j[n-2] = -1
j[n-1] = 1
j = (n-1) * j

lmbda = zeros(n)
L = transpose(D - B + A)
lmbda[1:n-1] = solve(L[1:n-1, 1:n-1], j[1:n-1])

#print "Lambda: ", lmbda
#print "L: \n", L

#print "L * lambda + j: ", dot(L, lmbda) + j

L[0,0] = +1.0
L[n-1,n-1] = +1.0
lmbda = solve(L, j)
print "Lambda: ", lmbda

fwd = D - B + A
fwd[0,0] = +1.0
fwd[n-1,n-1] = +1.0

dFda = zeros(n); dFda[0] = -1.0
duda = zeros(n)
for i in range(n):
  duda[i] = 1.0 - float(i) / (n-1)

print "dF/da: ", dFda
print "du/da: ", duda
print "fwd * du/da: ", dot(fwd, duda)
print "dF/da == -1 * fwd * du/da: ", dFda == -1 * dot(fwd, duda)

print "j * du/da: ", dot(j, duda)

print "lmbda * (fwd * du/da): ", dot(lmbda, dot(fwd, duda))
