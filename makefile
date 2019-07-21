CC = mpiCC
RM = rm
EXEC = MSR.a
LINKER = mpiCC
CFlags = -c -g -std=c++11
LFlags = -g
DEP  = petsc_real

# Edit include + lib dirs below
INCL = /usr/lib/petscdir/petsc3.10/x86_64-linux-gnu-real/include
LIB = /usr/lib/petscdir/petsc3.10/x86_64-linux-gnu-real/lib

install: Main.o Solver.o
	${LINKER} ${LFlags} Main.o Solver.o -l${DEP} -L${LIB} -o${EXEC}

Main.o: Main.cpp
	${CC} ${CFlags} Main.cpp -I${INCL}

Solver.o: Solver.cpp
	${CC} ${CFlags} Solver.cpp -I${INCL}

clean: Main.o Solver.o
	${RM} Main.o Solver.o ${EXEC}
