import os

import matplotlib.pyplot as plt
import numpy as np
import sympy as sp
from fluidity_tools import stat_parser

meshtemplate = """
Point(1) = {0, 0, 0, <dx>};
Extrude {0, 1, 0} {
  Point{1};Layers{<layers>};
}
Point(3) = {1, 0, 0, <dx>};
Extrude {0, 1, 0} {
  Point{3};Layers{<layers>};
}
Line(3)={1,3};
Line(4)={2,4};
Line Loop(5) = {4, -2, -3, 1};
Plane Surface(6) = {5};

Physical Line(1) = {1};
Physical Line(2) = {2};
Physical Line(3) = {4, 3};
Physical Surface(1) = {6};
"""


def generate_meshfile(name, layers):
    with open(name + ".geo", "w") as geo_file:
        geo_file.write(
            meshtemplate.replace("<dx>", str(1.0 / layers)).replace(
                "<layers>", str(layers)
            )
        )

    os.system("gmsh -2 " + name + ".geo")


def forcing(X):
    """Forcing function. Must be an analytic function of X[1] only"""

    return (X[1] ** 3, 0)


# Viscosity
mu = 1.0

# Note that because Coriolis can't be set from Python, the user has to ensure
# that this matches what it in the flml.
coriolis = 1.0


def analytic_solution(forcing):
    """Solve the ode d^2u/dx^2 = F/mu subject to u(0)=0, u(1)=0"""

    x = sp.Symbol("x")
    # Constants of integration.
    c1 = sp.Symbol("c_1")
    c2 = sp.Symbol("c_2")

    general = sp.integrate(sp.integrate(-forcing((0, x))[0] / mu, x) + c1, x) + c2

    constants = sp.solve(
        (sp.Eq(general.subs(x, 0), 0), sp.Eq(general.subs(x, 1), 0)), c1, c2
    )

    specific = general.subs(constants)

    return specific


def solution(forcing):
    """Return a function which is the solution to:
    d^2u/dx^2 = F/mu subject to u(0)=0, u(1)=0"""

    def sol(sx):
        return analytic_solution(forcing).subs(sp.Symbol("x"), sx[1])

    return sol


def analytic_pressure_solution(forcing):
    u = analytic_solution(forcing)

    return sp.integrate(-coriolis * u + forcing((0, sp.Symbol("x")))[1], sp.Symbol("x"))


def pressure_solution(forcing):
    """Return a function which is the solution to:
    dp/dx = f x u The constant of integration is set to 0."""

    def sol(sx):
        return analytic_pressure_solution(forcing).subs(sp.Symbol("x"), sx[1])

    return sol


def plot_theory():
    """Produce a plot showing the forcing, analytic velocity solution and
    analytic pressure solution"""

    plt.figure()

    y = np.linspace(0, 1, 21)

    psol = pressure_solution(forcing)

    usol = solution(forcing)

    v = 0 * y

    x = 0 * y

    us = np.array([float(usol(pos)) for pos in zip(x, y)])

    ps = np.array([float(psol(pos)) for pos in zip(x, y)])

    uf = np.array([forcing(pos) for pos in zip(x, y)])[:, 0]

    plt.subplots_adjust(wspace=0.25)
    plt.subplot(1, 3, 1)

    plt.quiver(x[1:-1], y[1:-1], uf[1:-1], v[1:-1], scale=1)
    plt.plot(uf, y)
    plt.xticks([0, 0.5, 1], map(str, [0, 0.5, 1]))
    plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1], map(str, [0, 0.2, 0.4, 0.6, 0.8, 1]))
    plt.ylabel("y")
    plt.xlabel("u source")

    plt.subplot(1, 3, 2)
    plt.plot(us, y)
    plt.quiver(x[1:-1], y[1:-1], us[1:-1], v[1:-1], scale=0.03)
    plt.xticks([0, 0.01, 0.02, 0.03], map(str, [0, 0.01, 0.02, 0.03]))
    plt.yticks([])
    plt.xlabel("u solution")

    plt.subplot(1, 3, 3)
    plt.plot(ps, y)
    plt.xticks([-0.02, -0.01, 0], map(str, [-0.02, -0.01, 0]))
    plt.yticks([])
    plt.xlabel("p solution")

    return uf, us, ps


def plot_results(dx, error):
    """plot_results(error)

    Produce a plot of the actual errors provided in the argument
    "error". Error should be a two column matrix with the first column being
    the velocity error and the second column the pressure error.
    """

    plt.figure()

    plt.loglog(dx, error)
    plt.loglog(dx, 0.03 * dx**2)
    plt.yticks(plt.yticks()[0], map(lambda x: "%3.1e" % x, plt.yticks()[0]))
    plt.xticks(plt.xticks()[0], map(lambda x: "%3.1e" % x, plt.xticks()[0]))

    plt.xlabel("dx")

    plt.title("Convergence of the rotating channel")

    plt.legend(("u error", "p error", "O(dx^2)"))


def retrieve_results(layers):
    """retrieve_results(layers)

    For each layer count in the sequence layers, retrieve the velocity and
    pressure error from the simulation results in appropriate channel-n
    directory.

    The first column of the result is the l2 norm of the error in the u
    component of velocity. The second is the l2 norm in the pressure.
    """
    from numpy import zeros

    error = zeros((len(layers), 2))

    for i, layer in enumerate(layers):
        s = stat_parser("channel-%d/rotating_channel.stat" % layer)

        error[i, 0] = s["Water"]["AnalyticUVelocitySolutionError"]["l2norm"][-1]
        error[i, 1] = s["Water"]["AnalyticPressureSolutionError"]["l2norm"][-1]

    return error
