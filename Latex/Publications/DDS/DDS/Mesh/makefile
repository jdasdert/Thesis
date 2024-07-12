FC = gfortran
FFLAGS = -llapack -g -fbounds-check 
SRC = createMesh.f90
SRC_old = createMesh_newOld.f90
MOD = tensorLib.f90 sortLib.f90

OBJ = $(SRC:.f90=.o)
OBJ_old = $(SRC_old:.f90=.o) 
MODOBJ = $(MOD:.f90=.o)
EXE = createMesh

DEBUGFLAGS = -g -llapack
DEBUGEXE = debugCreateMesh

mesh: $(MODOBJ) $(OBJ)
	$(FC) $(MODOBJ) $(OBJ) $(FFLAGS) -o $(EXE)

mesh_old: $(MODOBJ) $(OBJ_old)
	$(FC) $(MODOBJ) $(OBJ_old) $(FFLAGS) -o $(EXE)

debug: $(SRC) $(MOD) 
	$(FC) $(SRC) $(MOD) $(DEBUGFLAGS) -o $(DEBUGEXE)

test: test.f90
	$(FC) test.f90 -g -O2 -o test
	@./test


$(OBJ_old): $(MODOBJ) $(SRC_old)
	$(FC) -c -g $*.f90 -o $*.o

$(OBJ): $(MODOBJ) $(SRC)
	$(FC) -c -g $*.f90 -o $*.o

$(MODOBJ): $(MOD)
	$(FC) -c -g $*.f90 -o $*.o

clean:
	@rm -f *.mod *.o rm $(EXE)




