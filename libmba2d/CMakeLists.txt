set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ../lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ../lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../bin)

enable_language(Fortran)

add_library(mba2d STATIC
  ani2.f
  auxSE.f
  auxSF.f
  auxSP.f
  calcrv.f
  clpsF.f
  delaunay.f
  dsort.f
  error.f
  forlibmba.f
  insrtP.f
  isnan.c
  isnan_f.f
  lintrp2D.f
  list_bsct.f
  loadM.f
  makM.f
  makQ_linear.f
  makS.f
  mba_analytic.f
  mba_nodal.f
  minim.f
  moveP.f
  nlnfnc.f
  prjcrv.f
  saveM.f
  smoothing.f
  splitE.f
  swapF.f
  time.f
  untangle.f
  update.f
  utils.f
  ZZ.f
)
