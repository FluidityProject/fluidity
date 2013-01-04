from fluidity_tools import stat_parser as s
from math import fabs

def errors(stat1, stat2):
    s1=s(stat1)
    s2=s(stat2)
    l2error1=l2error(s1)
    l2error2=l2error(s2)
    inferr1=inferror(s1)
    inferr2=inferror(s2)
    return [l2error2/l2error1,inferr2/inferr1]

def l2error(stat):
    l2error=stat['Fluid']['LayerThicknessError']['l2norm'][-1]
    l2error=l2error/stat['Fluid']['LayerThickness']['l2norm'][0]
    return l2error

def inferror(stat):
    inferror=stat['Fluid']['LayerThicknessError']['max'][-1]
    inferror=max(inferror,fabs(stat['Fluid']['LayerThicknessError']['min'][-1]))
    inferror=inferror/stat['Fluid']['LayerThickness']['max'][0]
    return inferror

