from math import *
import array


def bottom(X,Y):
	b=[]
	for n in range(len(X)):
		pi=atan(1)*4
		b.append(400-100*exp(-pow((X[n]-0.5*max(X))/50000,2)))
#		b.append(400)
	return b


def U(X,Y):
	D=height(X,Y)
	u=[]
	for n in range(len(X)):
		u.append([])
		for m in range(2):
			u[n].append(0/D[n][m])	
	return u


def V(X,Y):
	v=[]
	for m in range(len(X)):
		v.append([])
		for n in range(2):
			v[m].append(0)	
	return v

def height(X,Y):
	D=[]
	b=bottom(X,Y)
	mx=max(X)
	for m in range(len(X)):
		D.append([])
		for n in range(2):
			if (n==0):
				D[m].append(100)
#				D[m][n] += \
#				   50*exp(-pow((X[m]-0.5*mx)/(0.1*mx),2))
			if (n==1):
				D[m].append(300)
				D[m][n] += b[m]-400
#				D[m][n] -= \
#				   50*exp(-pow((X[m]-0.5*mx)/(0.1*mx),2))
	return D


def middle(p1,p2):
	plist=[p1,p2]
	if (1 in plist) and (2 in plist):
		pr=6
	if (1 in plist) and (3 in plist):
		pr=5	
	if (3 in plist) and (2 in plist):
		pr=4
	return pr
