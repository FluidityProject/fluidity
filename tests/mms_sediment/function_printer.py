from mms_sediment_tools import *
from numpy import *
import matplotlib
import matplotlib.pyplot as plt
import sys

'''
run using:
python3 function_printer.py AA BB CC DD .. n_rows
where:
	AA, BB, CC, DD are names of functions in mms_rans_p2p1_keps_tools.py (any number can be entered)
	n_rows is the number of rows to display the functions on
'''

functions = []
for arg in sys.argv[1:-1]:
    functions.append(arg)
n_rows = int(sys.argv[-1])

matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'
plt.subplots_adjust(left=None, bottom=None, right=None, top=None,
                    wspace=0.2, hspace=0.2)
   
res = 50
X = linspace(0.0, pi, res)
Y = linspace(0.1, pi, res)
x = [0,0]

data = empty([len(functions), res, res])
for z, function in enumerate(functions):
    for j, x[0] in enumerate(X):
        for i, x[1] in enumerate(Y):
            data[z,i,j] = eval(function + '(x)')
    plt.subplot(n_rows, len(functions)/n_rows + 1, z+1)
    CS = plt.contour(X, Y, data[z])
    plt.clabel(CS, inline=1, fontsize=10, fmt='%1.1e')
    plt.title(functions[z])
    
plt.show()
