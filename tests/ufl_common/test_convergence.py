from fluidity_tools import stat_parser as stat
from math import log

def convergence_gal_stat(statA, statB):
    a_error = stat(statA)["Fluid"]["AbsoluteDifference"]["l2norm"][-1]
    b_error = stat(statB)["Fluid"]["AbsoluteDifference"]["l2norm"][-1]

    a_error_int = stat(statA)["Fluid"]["AbsoluteDifference"]["integral"][-1]
    b_error_int = stat(statB)["Fluid"]["AbsoluteDifference"]["integral"][-1]

    print a_error
    print b_error

    print a_error_int
    print b_error_int

    ab_ratio = a_error / b_error
    print "ratio", statA+"/"+statB+':', ab_ratio
    ab_ratio_int = a_error_int / b_error_int
    return [log(ab_ratio, 2), log(ab_ratio_int, 2)]

# vim:sw=4:ts=4:sts=4:et
