#!/usr/bin/env python3

import numpy
from numpy import array,argsort,corrcoef,size
import sys

def main():

    x = [0, 1, 2, 3, 4, 5, 6, 6.00000001]
    y = [0, 1, 2, 3, 4, 5, 6, 6.00000001]

    assert(interpolate(x, y, 5.5) == 5.5)
    assert(interpolate(x,y,0) ==0)
    assert(interpolate(x,y,6) ==6)
    assert(interpolate(x,y,4) ==4)  
    assert(interpolate(x,y,6.00000000001) == 6.00000000001)
    print("Success")

def interpolate(x, y, x0):
    """
    Interpolates array Y (y=f(x)) at point x0, returning y0
    x must be in increasing order
    """

    x_a = 0
    i_a = -1
    x_b = 0
    i_b = -1
    i = 0
    # do the literal edge cases - keep edges as the bounds n the y array
    if (x0 <= x[0]):
        return y[0]
    if (x0 >= x[-1]):
        return y[-1]

    for point in x:
        # out point exactly aligns with a point, simply return right value
        if (point == x0):
            return y[i]
        if (point < x0):
            x_a = point
            i_a = i
        if (point > x0):
            x_b = point
            i_b = i
            break
        i = i+1

    if (i_a == -1 or i_b == -1):
        sys.exit("Error interoplating")

    y0 = y[i_a] + (y[i_b]-y[i_a])*(x0-x_a)/(x_b-x_a)

    return y0


if __name__ == "__main__":
    main()
