from fluidity_tools import stat_parser as stat
from math import log

meshes = [['A','B'], ['B','C'], ['C','D']]

convergence_vel = [0,0,0]

print '-------------------------------------------------------------------------'
print 'Velocity convergence: l2Norm; integral; maximum'

for i, mesh in enumerate(meshes):

    a_error_x = stat("MMS_"+str(mesh[0])+".stat")["NS"]["Velocity_Difference%1"]["l2norm"][-1]
    b_error_x = stat("MMS_"+str(mesh[1])+".stat")["NS"]["Velocity_Difference%1"]["l2norm"][-1]
    a_error_y = stat("MMS_"+str(mesh[0])+".stat")["NS"]["Velocity_Difference%2"]["l2norm"][-1]
    b_error_y = stat("MMS_"+str(mesh[1])+".stat")["NS"]["Velocity_Difference%2"]["l2norm"][-1]
    
    a_error_int_x = stat("MMS_"+str(mesh[0])+".stat")["NS"]["Velocity_Difference%1"]["integral"][-1]
    b_error_int_x = stat("MMS_"+str(mesh[1])+".stat")["NS"]["Velocity_Difference%1"]["integral"][-1]
    a_error_int_y = stat("MMS_"+str(mesh[0])+".stat")["NS"]["Velocity_Difference%2"]["integral"][-1]
    b_error_int_y = stat("MMS_"+str(mesh[1])+".stat")["NS"]["Velocity_Difference%2"]["integral"][-1]
    
    a_error_inf_x = max([abs(stat("MMS_"+str(mesh[0])+".stat")["NS"]
                             ["Velocity_Difference%1"]["max"][-1]),
                         abs(stat("MMS_"+str(mesh[0])+".stat")["NS"]
                             ["Velocity_Difference%1"]["min"][-1])])
    b_error_inf_x = max([abs(stat("MMS_"+str(mesh[1])+".stat")["NS"]
                             ["Velocity_Difference%1"]["max"][-1]),
                         abs(stat("MMS_"+str(mesh[1])+".stat")["NS"]
                             ["Velocity_Difference%1"]["min"][-1])])
    a_error_inf_y = max([abs(stat("MMS_"+str(mesh[0])+".stat")["NS"]
                             ["Velocity_Difference%2"]["max"][-1]),
                         abs(stat("MMS_"+str(mesh[0])+".stat")["NS"]
                             ["Velocity_Difference%2"]["min"][-1])])
    b_error_inf_y = max([abs(stat("MMS_"+str(mesh[1])+".stat")["NS"]
                             ["Velocity_Difference%2"]["max"][-1]),
                         abs(stat("MMS_"+str(mesh[1])+".stat")["NS"]
                             ["Velocity_Difference%2"]["min"][-1])])

    ratio_x = a_error_x / b_error_x
    ratio_y = a_error_y / b_error_y
    ratio_int_x = a_error_int_x / b_error_int_x
    ratio_int_y = a_error_int_y / b_error_int_y
    ratio_inf_x = a_error_inf_x / b_error_inf_x
    ratio_inf_y = a_error_inf_y / b_error_inf_y

    convergence_vel[i] = [[log(ratio_x, 2), log(ratio_int_x, 2), log(ratio_inf_x, 2)], 
                          [log(ratio_y, 2), log(ratio_int_y, 2), log(ratio_inf_y, 2)]]

    print convergence_vel[i][0]
    print convergence_vel[i][1]

print '-------------------------------------------------------------------------'

fields = ['Pressure_Difference','Temperature_Difference','TurbulentKineticEnergy_Difference',
          'TurbulentDissipation_Difference']

convergence = [[0,0,0],[0,0,0],[0,0,0],[0,0,0]]

for i, field in enumerate(fields):
    print field + ' convergence: l2Norm; integral; maximum'

    for j, mesh in enumerate(meshes):
        
        a_error = stat("MMS_"+str(mesh[0])+".stat")["NS"][field]["l2norm"][-1]
        b_error = stat("MMS_"+str(mesh[1])+".stat")["NS"][field]["l2norm"][-1]
        
        a_error_int = stat("MMS_"+str(mesh[0])+".stat")["NS"][field]["integral"][-1]
        b_error_int = stat("MMS_"+str(mesh[1])+".stat")["NS"][field]["integral"][-1]
        
        a_error_inf = max([abs(stat("MMS_"+str(mesh[0])+".stat")["NS"]
                                 [field]["max"][-1]),
                             abs(stat("MMS_"+str(mesh[0])+".stat")["NS"]
                                 [field]["min"][-1])])
        b_error_inf = max([abs(stat("MMS_"+str(mesh[1])+".stat")["NS"]
                                 [field]["max"][-1]),
                             abs(stat("MMS_"+str(mesh[1])+".stat")["NS"]
                                 [field]["min"][-1])])
        
        ratio = a_error / b_error
        ratio_int = a_error_int / b_error_int
        ratio_inf = a_error_inf / b_error_inf
        
        convergence[i][j] = [log(ratio, 2), log(ratio_int, 2), log(ratio_inf, 2)]
        
        print convergence[i][j]
    
    print '-------------------------------------------------------------------------'
