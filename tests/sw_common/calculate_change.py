from vtktools import VtuDiff

def change(vtu1, vtu2):
    diff=VtuDiff(vtu1,vtu2)
    dh=diff.GetField('LayerThickness')
    maxdh=max(abs(dh.min()),abs(dh.max()))
    du=diff.GetField('Velocity')
    maxdu=max(abs(du.min()),abs(du.max()))
    return [maxdh, maxdu]
