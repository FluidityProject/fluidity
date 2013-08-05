import sys
import fluidity_tools
import vtktools

def norm_from_stat(filename,index_time_for_norm):
    st=fluidity_tools.stat_parser(filename)
    return st["ElapsedTime"]["value"][index_time_for_norm], st["fluid"]["ScalarAbsoluteDifference"]["l2norm"][index_time_for_norm], st["fluid"]["ScalarAbsoluteDifference"]["max"][index_time_for_norm]

filename = 'popbal_nonhomog_2d.stat'
timestep=0.005
outfile = open("time_versus_norm", "w")

for t in range(20):
	val=norm_from_stat(filename,((1.0/timestep)*(t+1))-1)
	outfile.write( "%g\t%g\t%g\n" %(val[0], val[1], val[2]) )

outfile.close()

