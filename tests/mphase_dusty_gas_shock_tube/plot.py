import matplotlib.pyplot as plt
import numpy
import shocktube
import vtktools

params = {
    "text.fontsize": 6,
    "legend.fontsize": 8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "lines.markersize": 6,
    "axes.titlesize": "small",
}
plt.rcParams.update(params)

# This corresponds to approximately t = 3.77e-4 seconds (normalised time of tau = 4 in
# the paper by Miura & Glass (1982))
filename = "mphase_dusty_gas_shock_tube_378.vtu"
vt = vtktools.vtu(filename)

# Time and coordinates
t = vt.GetScalarField("Air::Time")[0]
xyz = vt.GetLocations()
x = xyz[:, 0]

# Solution fields
p = vt.GetScalarField("Air::Pressure")
uvw_air = vt.GetVectorField("Air::Velocity")
u_air = uvw_air[:, 0]
uvw_dust = vt.GetVectorField("Dust::Velocity")
u_dust = uvw_dust[:, 0]
rho_air = vt.GetScalarField("Air::Density")
rho_dust = vt.GetScalarField("Dust::Density")
ie_air = vt.GetScalarField("Air::InternalEnergy")
ie_dust = vt.GetScalarField("Dust::InternalEnergy")
vfrac = vt.GetScalarField("Dust::PhaseVolumeFraction")

# Plot the numerical results.
# Note that we have divided the results by the reference
# velocity/pressure/density/internal energy to normalise them.
# x is normalised by the reference length.
plt.figure(1)
plt.subplot(321)
plt.plot(x / 0.0271002710027, p / 1.013e5, "b.", label="Numerical")

plt.subplot(322)
print(len(x))
print(len(u_air))
plt.plot(
    x[: len(x)] / 0.0271002710027, u_air / 286.980353992, "b.", label="Numerical (Air)"
)
plt.plot(
    x[: len(x)] / 0.0271002710027,
    u_dust / 286.980353992,
    "r.",
    label="Numerical (Particles)",
)

plt.subplot(323)
plt.plot(x / 0.0271002710027, rho_air / 1.23, "b.", label="Numerical")

plt.subplot(324)
plt.plot(x / 0.0271002710027, rho_dust, "b.", label="Numerical")
plt.title("Density of Particles")
plt.legend(loc=2)

plt.subplot(325)
plt.plot(x / 0.0271002710027, ie_air / 82357.7235772, "b.", label="Numerical (Air)")
plt.plot(x / 0.0271002710027, ie_dust / 82357.7235772, "r.", label="Numerical (Dust)")

# Multiply the PhaseVolumeFraction by the Dust Density to get the mass concentration,
# and then normalise.
plt.subplot(326)
plt.plot(x / 0.0271002710027, vfrac * 2500.0 / 1.23, "b.", label="Numerical")
plt.title("Mass Concentration from PhaseVolumeFraction of Particles")
plt.legend(loc=2)

# Now work out the frozen flow ('analytical') solutions
sol = numpy.array([shocktube.solution(xi, 4.0) for xi in x / 0.0271002710027])
p = sol[:, 0]
u = sol[:, 1]
rho = sol[:, 2]
ie = p / rho / (shocktube.gamma - 1.0)

plt.subplot(321)
plt.plot(x / 0.0271002710027, p, "-g", label="Frozen flow")
plt.title("Normalised Pressure")
plt.legend(loc=2)

plt.subplot(322)
plt.plot(x / 0.0271002710027, u, "g-", label="Frozen flow (Air)")
plt.title("Normalised Velocity")
plt.legend(loc=2)

plt.subplot(323)
plt.plot(x / 0.0271002710027, rho, "g-", label="Frozen flow")
plt.title("Normalised Density of Air")
plt.legend(loc=2)

plt.subplot(325)
plt.plot(x / 0.0271002710027, ie, "g-", label="Frozen flow")
plt.title("Normalised Internal Energy of Air")
plt.legend(loc=2)

plt.savefig("mphase_dusty_gas_shock_tube.png")

plt.show()
