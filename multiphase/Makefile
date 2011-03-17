
F90C = gfortran
LINKER = gfortran

#library options
AR = ar
ARFLAGS = cr

CC      = g++
CFLAGS	=

CXX      = g++
CXXFLAGS = 

# for normal compiling use
#F90FLAGS  = -I./include/  -fdefault-real-8 -O3  -ffast-math  -funroll-loops  -frecord-marker=4  -ffree-line-length-none  -ffixed-line-length-none 

# for compiling with debugging use
F90FLAGS  = -I./include/  -ffpe-trap=invalid,zero,underflow,overflow  -fbacktrace  -fdefault-real-8  -finit-integer=-10000  -finit-real=nan  -fbounds-check  -Wall -Wimplicit-interface -Wno-surprising -g -O0  -frecord-marker=4  -ffree-line-length-none  -ffixed-line-length-none 
                  
PROG_EXEC  =  ./bin/MultiphasePrototype.bin

PROG_LIB  =  ./lib/libMultiphasePrototype.a

LIBS = ${PROG_LIB} -ff2c -lc -lstdc++ -lgfortran -L/usr/local/lib -L/usr/lib 

PROG_OBJS = $(PROG_LIB)(multi_eos.o) \
            $(PROG_LIB)(multi_sparsity.o)  \
            $(PROG_LIB)(multi_matrix_operations.o) \
            $(PROG_LIB)(multi_shape_fct.o) \
            $(PROG_LIB)(cv-adv-dif.o) \
            $(PROG_LIB)(multi_dyncore_dg.o) \
            $(PROG_LIB)(multi_solvers.o) \
            $(PROG_LIB)(multi_compositional.o) \
            $(PROG_LIB)(multi_IO.o) \
            $(PROG_LIB)(Multiphase_Field_Advection.o) \
            $(PROG_LIB)(Multiphase_Mom_Press_Volf.o) \
            $(PROG_LIB)(Multiphase_Prototype.o)


VPATH=./source

.SUFFIXES: .cpp .F90 .o .a $(.SUFFIXES)
	                
.cpp.a:
	$(CXX) $(CXXFLAGS) -c $<
	$(AR) $(ARFLAGS) $@ $(<F:.cpp=.o)
	@rm -f $(<F:.cpp=.o)

.cpp.o:
	@echo "        CXX $<"
	$(CXX) $(CXXFLAGS) -c $<

.F90.a:
	$(F90C) $(F90FLAGS) -c $<
	$(AR) $(ARFLAGS) $@ $(<F:.F90=.o)

.F.a:
	$(F90C) $(F90FLAGS) -c $<
	$(AR) $(ARFLAGS) $@ $(<F:.F=.o)

# cancel implicit rule that tries to make .o from .mod with m2c:
%.o: %.mod       


default: 
	@mkdir -p ./lib ./include ./bin
	$(MAKE) $(PROG_LIB)
	$(MAKE) $(PROG_EXEC)

$(PROG_EXEC): $(PROG_LIB) main.o 
	ar x $(PROG_LIB) multi_prototype.o
	@touch dummy_for_mv_below.o
	@touch dummy_for_mv_below.mod
	@mv *.o ./include/.
	@mv *.mod ./include/.        
	$(LINKER) $(CXXFLAGS) -o $(PROG_EXEC) ./include/main.o $(LIBS) 

# clean everything, make code, make all test 
all:
	$(MAKE) superclean
	$(MAKE) default
	$(MAKE) test

# clean everything, make code, make all shorttest only 
allshort:
	$(MAKE) superclean
	$(MAKE) default
	$(MAKE) shorttest        


$(PROG_LIB):  $(PROG_OBJS)
$(PROG_LIB)(multi_eos.o): $(PROG_LIB)(multi_matrix_operations.o)
$(PROG_LIB)(multi_sparsity.o):
$(PROG_LIB)(multi_matrix_operations.o):
$(PROG_LIB)(multi_shape_fct.o):
$(PROG_LIB)(cv-adv-dif.o): $(PROG_LIB)(multi_solvers.o) $(PROG_LIB)(multi_shape_fct.o) $(PROG_LIB)(multi_matrix_operations.o) $(PROG_LIB)(multi_IO.o)
$(PROG_LIB)(multi_dyncore_dg.o): $(PROG_LIB)(multi_matrix_operations.o) $(PROG_LIB)(multi_shape_fct.o) $(PROG_LIB)(multi_solvers.o) $(PROG_LIB)(cv-adv-dif.o) $(PROG_LIB)(multi_IO.o)
$(PROG_LIB)(multi_solvers.o): $(PROG_LIB)(multi_shape_fct.o) $(PROG_LIB)(multi_matrix_operations.o)
$(PROG_LIB)(multi_compositional.o):
$(PROG_LIB)(multi_IO.o):
$(PROG_LIB)(Multiphase_Field_Advection.o): $(PROG_LIB)(multi_dyncore_dg.o) $(PROG_LIB)(multi_IO.o) $(PROG_LIB)(multi_compositional.o) 
$(PROG_LIB)(Multiphase_Mom_Press_Volf.o): $(PROG_LIB)(multi_dyncore_dg.o) $(PROG_LIB)(multi_eos.o) $(PROG_LIB)(multi_IO.o) 
$(PROG_LIB)(Multiphase_Prototype.o):	$(PROG_LIB)(multi_IO.o) $(PROG_LIB)(multi_sparsity.o) $(PROG_LIB)(Multiphase_Mom_Press_Volf.o) $(PROG_LIB)(Multiphase_Field_Advection.o)

# clean library/binary
clean:
	@rm -f $(PROG_EXEC) *.o *.mod ./lib/*.a

# clean library/binary and .o  
veryclean:
	@rm -f $(PROG_EXEC)  *.o *.mod ./include/*.o ./lib/*.a

# clean library/binary and .o .mod and tests   
superclean:
	@rm -f $(PROG_EXEC)  *.o *.mod  ./include/*.mod  ./include/*.o ./lib/*.a
	@cd tests; make clean

# commands for running test cases
# all the tests
test: tests

tests: shorttest mediumtest

# all the short tests
shorttest: shorttests

# all the short tests
mediumtest: mediumtests 

# the short tests - less than 2 min
shorttests:
	@cd tests; ../tools/testharness_MultiphasePrototype.py -l short

# the short tests - greater than 2 min, less than 10 min
mediumtests:
	@cd tests; ../tools/testharness_MultiphasePrototype.py -l medium

