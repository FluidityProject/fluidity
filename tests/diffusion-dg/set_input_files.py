import libspud
import shutil

class flml_things:
    "list of things to change in an flml file"

    def __init__(self,name='',mesh='',dt=0.2):
        self.name = name
        self.dt = dt
        self.mesh = mesh

parent_file = 'diffusion_dg1'
diffusion_dg_A = flml_things(name='diffusion_dg_A',mesh='MMS_A',dt=0.2)
diffusion_dg_B = flml_things(name='diffusion_dg_B',mesh='MMS_B',dt=0.2)
diffusion_dg_C = flml_things(name='diffusion_dg_C',mesh='MMS_C',dt=0.2)
diffusion_dg_D = flml_things(name='diffusion_dg_D',mesh='MMS_D',dt=0.2)
diffusion_dg_E = flml_things(name='diffusion_dg_E',mesh='MMS_E',dt=0.2)
filelist = [diffusion_dg_A,diffusion_dg_B,diffusion_dg_C,diffusion_dg_D,diffusion_dg_E]

for file in filelist:
    shutil.copy(parent_file+'.flml',file.name+'.flml')
    libspud.load_options(file.name+'.flml')
    libspud.set_option('/simulation_name', file.name)
    libspud.set_option('/geometry/mesh[@name="CoordinateMesh"]/from_file/file_name', file.mesh)
    libspud.set_option('/timestepping/timestep', file.dt)

    libspud.write_options(file.name+'.flml')
