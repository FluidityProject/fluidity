#!/usr/bin/env python

import fluidity.diagnostics.annulus_mesh as mesh
import fluidity.diagnostics.triangletools as tt

div = mesh.SliceCoordsConstant(0.0, 1.0, 3)
m = mesh.GenerateRectangleMesh(div, div)
tt.WriteTriangle(m, "square-structured-linear")
