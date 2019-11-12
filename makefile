CC = mpiCC
RM = rm
EXEC = MSR.a
LINKER = mpiCC
CFlags = -c -g -std=c++11
LFlags = -g
DEP  = petsc_real

# Edit include + lib dirs below
INCL = /usr/lib/petscdir/3.10/include
LIB = /usr/lib/petscdir/3.10/lib

install: main.o solver.o inout.o
	${LINKER} ${LFlags} main.o solver.o inout.o -l${DEP} -L${LIB} -o${EXEC}

main.o: src/main.cpp
	${CC} ${CFlags} src/main.cpp -I${INCL}

solver.o: src/solver.cpp
	${CC} ${CFlags} src/solver.cpp -I${INCL}

inout.o: src/inout.cpp
	${CC} ${CFlags} src/inout.cpp -I${INCL}

clean: main.o solver.o inout.o
	${RM} main.o solver.o inout.o ${EXEC}
