import vtk
import glob
import operator
import numpy
import pickle
from scipy.interpolate import interp1d

def GetFileList(ele_type):
    if ele_type[-2:]=='DG':
        CG=False
    else:
        CG=True

    filelist=glob.glob('Sloshing_Water_Tank_%s_*.vtu'%(ele_type))
    nums=[]
    for i in range(len(filelist)):
        nums.append([i, int(filelist[i].split(".vtu")[0].split("_")[-1])])
    nums.sort(key=operator.itemgetter(1))
    return [filelist[N] for N in numpy.array(nums)[:,0]]

def GetRunData(ele_type):
    data=[]
    for k,file in enumerate(GetFileList(ele_type)):
        X,Y=GetContour(file)
        data.append((X,Y))
    return data

def SetStoredData(ele_types):
    for ele_type in ele_types:
        file=open('%s_regression.data'%ele_type,'w')
        data=GetRunData(ele_type)
        pickle.dump(data,file)
        file.close()


def TestContour(ele_types):
    flags=[]
    stored_types=GetStoredData(ele_types)
    for ele_type,stored_type in zip(ele_types,stored_types):
        flag=True
        data=GetRunData(ele_type)
        for d,s in zip(data,stored_type):
            flag= flag and LInfContourDifference(d,s)<1.0e-2
        flags.append(flag)
    return flags

def GetStoredData(ele_types):
    types=[]
    for ele_type in ele_types:
        file=open('%s_regression.data'%ele_type,'r')
        data=pickle.load(file)
        types.append(data)
        file.close()
    return types

def GetContour(file,val=0.5):
    reader=vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(file)
    reader.Update()
    data=reader.GetOutput()
    data.GetPointData().SetActiveScalars("Component1::ComponentMassFractionPhase1")
    contour=vtk.vtkContourFilter()
    contour.SetInput(data)
    contour.SetValue(0,val)
    contour.Update()
    X=[]
    Y=[]
    polydata=contour.GetOutput()
    for n in range(polydata.GetNumberOfPoints()):
        p=polydata.GetPoints().GetPoint(n)
        X.append(p[0])
        Y.append(p[1])
    X=numpy.sort(numpy.array(X))
    Y=numpy.sort(numpy.array(Y))
    return X,Y

def ContourDifference(C1,C2):
    X=numpy.concatenate((C1[0],C2[0]))
    X=numpy.sort(X)
    c1=interp1d(C1[0],C1[1])
    c2=interp1d(C2[0],C2[1])

    y=c1(X)-c2(X)
    return X,y


def LInfContourDifference(C1,C2):
    X,Y=ContourDifference(C1,C2)
    return abs(Y).max() 


    
    
