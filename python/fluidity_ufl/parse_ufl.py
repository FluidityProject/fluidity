        
class Function():
    def __init__(self, field):

        self.rank=field.rank
        self.test=False
        self.trial=False

        self.field=field
        self.name=field.name

        self.rank=field.rank
        self.indices=[slice(None,None,None) for i in range(field.rank)]
        self.form_gradient=False

    def __mul__(self, other):

        if isinstance(other, Form):
            return Form.__rmul__(self)
        else:
            product=Form()
            product=product.__mul__(self)
            product=product.__mul__(other)
            return product

    def __getitem__(self, i):

        try:
            if (len(i)!=self.rank):
                raise TypeError("Wrong number of subscripts")

            for ii in i:
                if (not isinstance(ii, slice)) and (not isinstance(ii, index)):
                    raise TypeError("Indices must be of type slice or index.")

            self.indices=list(i)

        except AttributeError:
            # In this case there is only one argument.
            if (self.rank!=1):
                raise TypeError("Wrong number of subscripts")

            if (not isinstance(i, slice)) and (not isinstance(i, index)):
                raise TypeError("Indices must be of type slice or index.")

            self.indices=[i]

        return self

    def current_rank(self):
        """Return the rank of a function taking into account any reduction
        of rank caused by repeated indices."""

        rank=self.rank

        indices={}
        for i in self.indices:
            if (isinstance(i, index)):
                if i.id in indices:
                    indices[i.id]=indices[i.id]+1

        for i in indices.itervalues():
            if i==2:
                rank=rank-2
            if i>2:
                raise TypeError("Index repeated more than twice")
        
class BasisFunction(Function):     
    pass

class TestFunction(BasisFunction):
    def __init__(self, field):
        Function.__init__(self,field)

        self.name=self.name+"Test"

        self.test=True

class TrialFunction(BasisFunction):
    def __init__(self, field):
        Function.__init__(self,field)

        self.name=self.name+"Trial"

        self.trial=True

class grad(Function):
    def __init__(self, primitive):
        Function.__init__(self,primitive.field)

        self.test=primitive.test
        self.trial=primitive.trial

        self.rank=primitive.rank+1
        self.indices=[slice(None,None,None)]+self.indices

        self.name="grad_"+primitive.name

        # Note the nasty side effect!
        primitive.form_gradient=True

        
        grad.primitive=primitive

class Measure():
    def __init__(self, name):
        self.name=name

    def __mul__(self, other):

        if isinstance(other, form):
            return form.__mul__(self)
        else:
            product=Form()
            product=product.__mul__(self)
            product=product.__mul__(other)
            return product

index_count=0

class index():
    def __init__(self):
        global index_count
        
        index_count=index_count+1

        self.id=index_count

i=index()
j=index()
k=index()

dx=Measure("dx")

state=state_type()

class Integral():
    def __init__(self, source=None ):

        self.test=None
        self.trial=None
        self.measure=None

        self.rank=0
        self.sum_indices={}
        self.functions=[]
        self.scalars=[]

        if isinstance(source, type(None)):
            pass
        elif isinstance(source, Integral):
            self=source.copy()
            
        elif isinstance(source, Function):
            if source.test:
                self.test=source
            if source.trial:
                self.trial=source

            self.functions.append(source)

            for ii in source.indices:
                if (isinstance(ii, index)):
                    if ii.id in self.sum_indices:
                        self.sum_indices[ii.id]= \
                                            product.sum_indices[ii.id]+1
                    else:
                        self.sum_indices[ii.id]=1
                elif (isinstance (ii, slice)):
                    self.rank=product.rank+1
            
            
        elif isinstance(source, Measure):
            self.measure=source

        else:
            raise TypeError("Unable to form Integral from "+str(source))
        
    def copy(self):
        copy=Integral()

        copy.test=self.test
        copy.trial=self.trial
        copy.measure=self.measure

        copy.rank=self.rank
        copy.sum_indices=self.sum_indices        
        copy.functions=self.functions
        copy.scalars=self.scalars

        return copy
        
    def invalid(self):
        '''Decide if self is a valid linear or bilinear form. If valid,
        return none. Otherwise return a string describing what has gone wrong.'''

        str=""

        if (not (self.rank in (0,1,2))):
            str=str+"Rank is %d but must be 0, 1 or 2.\n"%self.rank

        if (not (self.test)):
            str=str+"There is no test function.\n"

        if (not (self.measure)):
            str=str+"There is no measure to integrate over.\n"

        for ii in self.sum_indices.iteritems():
            if ii[1]<2:
                str=str+"Index %d only occurs once.\n"

        if (str==""):
            return None
        else:
            return str
        
    def __neg__(self):
        product=self.copy()
        product.scalars.append(-1)

    def __rmul__(self, other):

        # Multiplication is commutative!
        return self.__mul__(other)
            
    def __mul__(self, other):
        product=self.copy()

        if (isinstance(other, Function)):
            if other.test:
                if product.test:
                    raise TypeError("An integral must have exactly one test function")
                else:
                    product.test=other
                    
            if other.trial:
                if product.trial:
                    raise TypeError("An integral may have at most one trial function")
                else:
                    product.trial=other

            product.functions.append(other)

            for ii in other.indices:
                if (isinstance(ii, index)):
                    if ii.id in product.sum_indices:
                        product.sum_indices[ii.id]= \
                                            product.sum_indices[ii.id]+1
                    else:
                        product.sum_indices[ii.id]=1
                elif (isinstance (ii, slice)):
                    product.rank=product.rank+1

        elif (isinstance(other, Measure)):
            if product.measure:
                raise TypeError("An integral must have exactly one measure")

            product.measure=other

        elif (isinstance(other, float) or isinstance(other, int)):
            product.scalars.append(other)

        else:
            raise TypeError("Cannot multiply "+str(self)+" * "+str(other))

        return product

class Form(list):
    """This is in effect just a list of integrals."""

    def __sub__(self, other):
        if not isinstance(other, Form):
            raise TypeError("Subtraction of forms is only supported between forms")

        return self+(-other)
        

    def __mul__(self, other):
        product=Form()

        if (len(self)==0):
            product.append(Integral(other))
        else:
            for i in self:
                product.append(i*other)

        return product

    def __rmul__(self, other):
        product=Form()

        if (len(self)==0):
            product.append(Integral(other))
        else:
            for i in self:
                product.append(other*i)

        return product

    def __add__(self,other):
        sum=Form()

        for s in self:
            sum.append(s)
        for o in other:
            sum.append(o)

        return sum
    

    def copy(self):
        copy=Form()
        for i in self:
            copy.append(i.copy)


class field():
    def __init__(self, name):
        self.name=name
        field.rank=0

        self.addtos=[]
        self.solve=None
        
    def __iadd__(self, other):

        if x.solve:
            raise AttributeError("Cannot assemble x by both addition and solve.")


        if not isinstance(other, Form):
            raise TypeError("Only know how to add forms to fields.")

        for integral in other:
            if (integral.trial):
                raise AttributeError("Cannot add a form with a trial function "+\
                                         "to a field.")
            if integral.rank!=field.rank:
                raise AttributeError("Rank mismatch in field addto")
            self.addtos.append(integral)


class scalarField(field):
    def __init__(self, name):
        field.__init__(self, name)
        
class vectorField(field):
    def __init__(self, name):
        field.__init__(self, name)
        self.rank=1
        
class tensorField(field):
    def __init__(self, name):
        field.__init__(self, name)
        self.rank=2

class __fake_field_dict__():
    def __init__(self,rank):
        self.rank=rank
        
    def __getitem__(self,name):
        if (self.rank==0):
            return scalarField(name)
        elif(self.rank==1):
            return vectorField(name)
        elif(self.rank==2):
            return tensorField(name)
        endif
        
        
class state_type():
    def __init__(self):
        self.scalarFields=__fake_field_dict__(0)
        self.vectorFields=__fake_field_dict__(1)
        self.tensorFields=__fake_field_dict__(2)

class Matrix():
    def __init__(self, test, trial):

        if not isinstance(test, TestFunction):
            raise TypeError("The test function of matrix must have type TestFunction")

        if not isinstance(trial, TrialFunction):
            raise TypeError("The trial function of matrix must have type TrialFunction")

        self.test=test
        self.trial=trial

        self.addtos=[]

    def __iadd__(self, other):

        if not isinstance(other, Form):
            raise TypeError("Only know how to add forms to matrices.")

        for integral in other:
            if not (integral.trial):
                raise AttributeError("Cannot add a form with no trial function "+\
                                         "to a matrix.")
            if (integral.test.field.name!=self.test.field.name):
                raise AttributeError("Test function mismatch.")
            if (integral.trial.field.name!=self.trial.field.name):
                raise AttributeError("Test function mismatch.")
            self.addtos.append(integral)

        return self

def solve(x, A, b):
    """Solve the matrix equation Ax=b where A is a matrix and b, x are
    fields."""

    if not isinstance(A, Matrix):
        raise TypeError("A must be a matrix.")
    if not isinstance(b, field):
        raise TypeError("b must be a field.")
    if not isinstance(x, field):
        raise TypeError("x must be a field.")
    
    if x.addtos:
        raise AttributeError("Cannot assemble x by both addition and solve.")
    if x.solve:
        raise AttributeError("Already have a solution algorithm for x.")

    x.solve=(A,b)
    
def parse_file(file):

    global_dict={}
    local_dict={}
    
    exec("from parse_ufl import *", global_dict)

    execfile(file, global_dict, local_dict)

    return(local_dict)

class formsum(list):
    pass
