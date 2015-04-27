#!/usr/bin/env python

import fluidity.diagnostics.annulus_mesh as mesh
import fluidity.diagnostics.gmshtools as gt

div = mesh.SliceCoordsConstant(0.0, 1.0, 3)
m = mesh.GenerateRectangleMesh(div, div)
gt.WriteMsh(m, "square-structured-linear.msh")
