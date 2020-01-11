Reading H5Part data into VisIt:

VisIt, http://www.llnl.gov/visit,  is an open source point-and-click 3D scientific visualization application that supports most of the common visualization techniques on structured and unstructured grids. One of its advantages is that it employs a distributed and parallel architecture in order to handle extremely large data sets interactively. VisIt's rendering and data processing capabilities are split into viewer and engine components that may be distributed across multiple machines.

This database plugins allows the user to read H5Part data into VisIt. To be recognized by VisIt the filename needs to have the extension ".h5part". This is important, otherwise VisIt will try to open the file with another reader and it will fail.

Building the plugin:
You must first install VisIt in your machine to be able to link to its libraries. 

Modify the Makefile to reflect your installation.

TOPDIR is the path to the VisIt Distribution installation.
HOMEPLUGINS is the path to your home dir .visit plugin directory. Create it if not there.
SHLIB_FORCED is the path to your H5Part. I made H5Part statically linked to the hdf5 library, if you don't, add your -L/hdf5librarypath -lhdf5 
Modify the location of the hdf5 include and the H5Part library.


If compiling a parallel version add -DPARALLEL_IO to the CXXFLAGS and CPPFLAGS.


Note: TOPDIR is defined twice because in make-variables it is defined by the VisIt developers to be in one of their machines.
