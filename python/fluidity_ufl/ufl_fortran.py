
import parse_ufl

detwei=("detwei", "detwei", "real, dimension(ngi), intent(in) :: detwei")

class fortran_module():
    def __init__(self, filename):
        """Given a ufl file produce a fortran module implementing that
        ufl. The module will have the name basename_module where basename is
        filename with any trailing .ufl stripped. The module will have a
        single public subroutine named basename. Basename will take a
        state_type as an argument."""

        self.filename=filename
        if filename[-4:]==".ufl":
            self.basename=filename[:-4]
        else:
            self.basename=filename

        self.parsed_ufl=parse_ufl.parse_file(filename)

        

class fortran_integral():
    def __init__(self, integral, dim="dim"):
        err=integral.invalid()
        if err:
            raise TypeError(err)

        self.integral=integral
        self.test_dim="ele_loc("+self.integral.test.name+")"
        try:
            self.trial_dim="ele_loc("+self.integral.trial.name+")"
        except AttributeError:
            self.trial_dim=None

        # Dimension of the problem space.
        self.dim=dim
        
        self.set_arguments()
        self.map_sum_indices()

    def fortran(self):

        dim_indices=["dim"+str(i)+"_i" for i in range(self.integral.rank)]

        name, declaration, lhs=self.function_spec(dim_indices)
        
        decs=self.argument_declarations+[declaration]

        core_indices="integer :: gi, test_i"
        if self.integral.trial:
            core_indices=core_indices+", trial_i"

        decs=decs+[core_indices]

        if dim_indices:
            decs=decs+["integer :: "+",".join(dim_indices)]

        if self.sum_index_count:
            sum_indices="integer :: "+\
                  ", ".join(["i"+str(i+1) for i in
                             range(self.sum_index_count)])
            decs=decs+[sum_indices]

        body_lhs=[lhs+" = "+lhs+"&"]
        body = "+ "
        # Deep copy of dim_indices to be consumed in the next loop.
        ldim_indices=list(dim_indices)
        for f in self.iterfunctions():
            body=body+self.index_function(f, ldim_indices)+"*"
        body=body+"detwei(gi)"
        body=body_lhs+indent([body])

        # Note that the innermost loops comes first.

        for index in dim_indices:
            body=do_loop(index,self.dim,body)
        for index in range(self.sum_index_count):
            body=do_loop("i"+str(index+1),self.dim,body)
        
        body=do_loop("gi","ngi",body)
        if self.integral.trial:
            body=do_loop("trial_i",self.trial_dim,body)
        body=do_loop("test_i",self.test_dim,body)
        
        
        code=decs
        code=code+[""]
        code=code+["integral=0.0"]
        code=code+[""]
        code=code+body
        code=code+[""]

        code=indent(code)
        
        code=["function "+name+"("+", ".join(self.dummy_arguments)+")"\
              +" result (integral)"]+code

        code=code+["end function "+name]

        return code

    def index_function(self,f, dim_indices):
        """ Produce a string containing the name of f and the appropriate
        indices"""

        index_code=[]
        for i in f.indices:
            if (isinstance(i, slice)):
                index=dim_indices[0]
                dim_indices.remove(index)
                index_code.append(index)
            else:
                index="i"+str(self.sum_index_map[i.id])
                index_code.append(index)

        # Now the special indices for basis functions.
        if f.trial:
            index_code.append("trial_i")
        if f.test:
            index_code.append("test_i")

        # We always need to integrate over quadrature points.
        index_code.append("gi")
        
        code=f.name+"("+",".join(index_code)+")"

        return code
        
    def iterfunctions(self):
        """Generator enabling iteration of the the functions of a
        fortran_integral in the order they appear in the argument list."""

        yield self.integral.test

        for f in self.integral.functions:
            # We've already done test and trial.
            if f.test or f.trial:
                continue

            yield f

        if self.integral.trial:
            yield self.integral.trial


    def function_spec(self, dim_indices):

        name="integral"
        declaration="real, dimension("
        lhs_args=dim_indices+["test_i"]

        for f in self.iterfunctions():
            name=name+"_"+f.name
            
            for i in f.indices:
                # Slices are free indices.
                if isinstance(i, slice):

                    declaration=declaration+str(self.dim)+", "
                    
                else:

                    name=name+"_i"+str(self.sum_index_map[i.id])

        name=name+"_"+self.integral.measure.name

        declaration=declaration+self.test_dim+", "

        if self.integral.trial:
            declaration=declaration+self.trial_dim+", "
            lhs_args=lhs_args+["trial_i"]

        declaration=declaration+"ngi) :: integral"
        lhs="integral("+", ".join(lhs_args)+")"
        
        return name, declaration, lhs

    def map_sum_indices(self):
        """Form a list of local index names corresponding to the global
        indices recorded in the integral"""

        self.sum_index_map={}
        self.sum_index_count=0

        for ii in self.integral.sum_indices.iterkeys():
            if ii not in self.sum_index_map:
                self.sum_index_count=self.sum_index_count+1
                
                self.sum_index_map[ii]=self.sum_index_count
        

    def set_arguments(self):
        self.argument_declarations=[]
        self.variable_declarations=[]
        
        self.dummy_arguments=[]
        self.actual_arguments=[]

        for f in self.iterfunctions():

            args=function_to_arguments(f)
            self.dummy_arguments.append(args[0])
            self.actual_arguments.append(args[1])
            self.argument_declarations.append(args[2])

        # You always need detwei.
        self.dummy_arguments.append(detwei[0])
        self.actual_arguments.append(detwei[1])
        self.argument_declarations.append(detwei[2])
    
def function_to_arguments(function):
    '''Take a function object and return actual and dummy argument names as
    well as a declaration of the dummy argument'''

    actual=function.name
    dummy=function.name

    if isinstance(function, parse_ufl.BasisFunction):
        declaration="element_type, intent(in) :: "+dummy
    else:
        if function.test or function.trial:
            extra_args=2
        else:
            extra_args=1
        declaration="real, dimension("+",".join([":" for i in
                                                 range(function.rank+extra_args)])\
                     +"), intent(in) :: "+dummy

    return actual, dummy, declaration

def indent(code):
    
    return ["  "+line for line in code]

def do_loop(var, size, body):

    code=["do "+var+" = 1, "+str(size)]\
          + indent(body)\
          +["end do"]
    
    return code
