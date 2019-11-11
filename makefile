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

install: Main.o Solver.o
	${LINKER} ${LFlags} Main.o Solver.o -l${DEP} -L${LIB} -o${EXEC}

Main.o: src/Main.cpp
	${CC} ${CFlags} src/Main.cpp -I${INCL}

Solver.o: src/Solver.cpp
	${CC} ${CFlags} src/Solver.cpp -I${INCL}

clean: src/Main.o src/Solver.o
	${RM} src/Main.o src/Solver.o ${EXEC}
