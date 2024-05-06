from math import log

from fluidity_tools import stat_parser as stat


def report_convergence(file1, file2):
    print(file1, "->", file2)

    stat1 = stat(file1)
    stat2 = stat(file2)

    errortop_l2_1 = stat1["Fluid"]["FreeSurfaceDifference"]["surface_l2norm%Top"][-1]
    errortop_l2_2 = stat2["Fluid"]["FreeSurfaceDifference"]["surface_l2norm%Top"][-1]
    convergencetop_l2 = log((errortop_l2_1 / errortop_l2_2), 2)

    print("  convergencetop_l2 = ", convergencetop_l2)
    print("    errortop_l2_1 = ", errortop_l2_1)
    print("    errortop_l2_2 = ", errortop_l2_2)

    errorbottom_l2_1 = stat1["Fluid"]["FreeSurfaceDifference"]["surface_l2norm%Bottom"][
        -1
    ]
    errorbottom_l2_2 = stat2["Fluid"]["FreeSurfaceDifference"]["surface_l2norm%Bottom"][
        -1
    ]
    convergencebottom_l2 = log((errorbottom_l2_1 / errorbottom_l2_2), 2)

    print("  convergencebottom_l2 = ", convergencebottom_l2)
    print("    errorbottom_l2_1 = ", errorbottom_l2_1)
    print("    errorbottom_l2_2 = ", errorbottom_l2_2)

    error_l2_1 = stat1["Fluid"]["FreeSurfaceDifference"]["surface_l2norm%Both"][-1]
    error_l2_2 = stat2["Fluid"]["FreeSurfaceDifference"]["surface_l2norm%Both"][-1]
    convergence_l2 = log((error_l2_1 / error_l2_2), 2)

    print("  convergence_l2 = ", convergence_l2)
    print("    error_l2_1 = ", error_l2_1)
    print("    error_l2_2 = ", error_l2_2)

    error_linf_1 = stat1["Fluid"]["FreeSurfaceDifference"]["max"][-1]
    error_linf_2 = stat2["Fluid"]["FreeSurfaceDifference"]["max"][-1]
    convergence_linf = log((error_linf_1 / error_linf_2), 2)

    print("  convergence_linf = ", convergence_linf)
    print("    error_linf_1 = ", error_linf_1)
    print("    error_linf_2 = ", error_linf_2)

    return [convergencetop_l2, convergencebottom_l2, convergence_linf]
