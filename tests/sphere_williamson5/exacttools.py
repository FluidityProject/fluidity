import numpy

def setexact(coords_filename,solution_filename,dim=3,tol=1.0e-8):
    #return down
    from numpy import fromfile, reshape, where
    Coords = fromfile(coords_filename)
    Coords = reshape(Coords,(Coords.size/dim,dim))
    Vals = fromfile(solution_filename)
    print Coords.shape, Vals.size
    def val(X,t):
        print X
        #Search amongst Coords for X
        B = ((Coords[:,0]-X[0])**2+
             (Coords[:,1]-X[1])**2+
             (Coords[:,2]-X[2])**2)**0.5
        vals = where(B<tol,Vals)
        assert(vals.size==1)
        return vals[0]
    return val

def makecoordsfile(filename):
    nodefile = file(filename+'.node', 'r')
    row1 = numpy.double(nodefile.readline().split())
    assert(row1.shape==(4,))
    nnodes=int(row1[0])
    dim=int(row1[1])
    assert(row1[2]==0)
    assert(row1[3]==0)
    X = numpy.zeros((nnodes,dim))
    for node in range(nnodes):
        row = numpy.double(nodefile.readline().split())
        assert(row.shape==(dim+1,))
        X[node,:] = row[1:dim+1]
    X.tofile(filename+'.dat',sep=' ')

def makedat(filename,dim=3):
    a = numpy.fromfile(filename+'.node',sep=' ')
    n = a.size
    a = numpy.reshape(a,(n/dim,dim))
    a.tofile(filename+'.dat',sep=' ')
    
