#!/usr/bin/env python

import numpy
from numpy import array,argsort,corrcoef,size
import sys

def main():

    x = [0, 1, 2, 3, 4, 5, 6]
    y = [0, 1, 2, 3, 4, 5, 6]

    print(interpolate(x, y, 5.5))
    assert(interpolate(x, y, 5.5) == 5.5)
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
    for point in x:
        if (point <= x0):
            x_a = point
            i_a = i
        if (point > x0):
            x_b = point
            i_b = i
            break
        i = i+1

    if (i_a == -1 or i_b == -1):
        sys.exit("Error interpolating")

    y0 = y[i_a] + (y[i_b]-y[i_a])*(x0-x_a)/(x_b-x_a)

    return y0


if __name__ == "__main__":
    main()
