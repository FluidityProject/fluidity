import vtktools


def values_per_node(file):
    u = vtktools.vtu(file)
    zoo = u.GetScalarField("Zooplankton")
    phyto = u.GetScalarField("Phytoplankton")
    nut = u.GetScalarField("Nutrient")
    det = u.GetScalarField("Detritus")

    return phyto, zoo, nut, det
