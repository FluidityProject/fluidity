IMPORTANT NOTES ON BACKWARD FACING STEP 3D EXAMPLE:

Please note, this example is intended to be run in parallel.
Before running this test, ensure that you have run 'make fltools' in
the fluidity directory in order to build the fldecomp tool. This tool is
needed to decompose the mesh into parallel domains, and is called
when you run 'make preprocess' in this directory.

The example is currently set to run on 4 processors, so if you wish to
run it on a different number (depending on your machine) then alter the
Makefile accordingly. The relevant lines are:
<../../bin/fldecomp -n 4 -f step3d>
where <-n [number]> specifies the number of mesh parts, and:
<mpiexec -np 4 ../../bin/fluidity -v2 -l backward_facing_step_3d.flml>
where <-np [number]> specifies the number of processors on which to run.
The numbers must be the same or the executable will not run.

When running Fluidity in parallel, the output file format automatically
changes to .pvtu, with each processor outputting .vtu files on its
allocated mesh fragment (these are saved in subdirectories).
Pvtu files open with Mayavi etc. exactly the same way as vtu files.
However the postprocess script only works for .pvtu files, so if you
want to run on 1 processor the script will not work.

For more details see entries for "parallel" in the index of the Fluidity manual,
and the 'Backward Facing Step' subsection in Chapter 10 (Examples).
