#!/usr/bin/env python3

import fluidity.diagnostics.annulus_mesh as annulus_mesh
import fluidity.diagnostics.triangletools as triangletools

coordsY = annulus_mesh.SliceCoordsLinear(0.0, 1.0, 0.01, 20)
coordsX = annulus_mesh.SliceCoordsLinear(0.0, 1.0, 0.01, 20)
mesh = annulus_mesh.GenerateRectangleMesh(coordsX, coordsY)
triangletools.WriteTriangle(mesh, "square-structured-linear")

