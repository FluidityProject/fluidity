import numpy

def make_node_file(filename,X,Y,Z):
    print "Creating .node file "+filename+"."

    f = open(filename+'.node','w')
    assert(numpy.shape(X)==numpy.shape(Y))
    assert(numpy.shape(X)==numpy.shape(Z))
    f.write(str(numpy.size(X))+' '+'3 0 0\n')

    for vert in range(numpy.size(X)):
        x = X[vert]
        y = Y[vert]
        z = Z[vert]
        f.write(str(vert+1)+' '+str(x)+' '+str(y)+' '+str(z)+'\n')

    f.write('Produced by: make_test_manifold_meshes')
    f.close()
    print ".node file Complete."

def make_ele_file(filename,EVList):
    print "Creating .ele file "+filename+"."
    elements = numpy.size(EVList,0)
    f = open(filename+'.ele','w')
    f.write(str(elements)+' 3 1\n')
    for ele in range(elements):
        f.write(str(ele+1)+' ')
        f.write(str(int(EVList[ele,0]))+' ')
        f.write(str(int(EVList[ele,1]))+' ')
        f.write(str(int(EVList[ele,2]))+' 0\n')
    f.write('Produced by: make_test_manifold_meshes')
    f.close()
    print ".ele file Complete."

print "Making flat mesh"

#numbering is  43
#              12
#and 5 in the middle.

X = numpy.array([0.0,1.0,1.0,0.0,0.5])
Y = numpy.array([0.0,0.0,1.0,1.0,0.5])
Z = numpy.array([0.,0.,0.,0.,0.])
EVList = numpy.zeros((4,3),dtype=numpy.integer)
EVList[0,:] = numpy.array([1,2,5])
EVList[1,:] = numpy.array([2,3,5])
EVList[2,:] = numpy.array([3,4,5])
EVList[3,:] = numpy.array([4,1,5])

make_node_file('flat',X,Y,Z)
make_ele_file('flat',EVList)

print "Making sloped mesh"

Z = X-0.5
make_node_file('sloped',X,Y,Z)
make_ele_file('sloped',EVList)

print "Making random height mesh"

Z = numpy.random.rand(numpy.size(X))-0.5
make_node_file('rand',X,Y,Z)
make_ele_file('rand',EVList)

