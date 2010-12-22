xyz2grd Bathymetry.xyz -GBathymetry.grd -R0/5.488/0/3.402 -I0.014
grdmath Bathymetry.grd -1.0 MUL = Bathymetry.grd
