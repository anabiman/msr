CC = mpiCC
RM = rm
EXEC = MSR.a
LINKER = mpiCC
CFlags = -c -g -std=c++0x
LFlags = -g
INCL1 = ${HOME}/.local/usr/include
INCL2 = ${HOME}/.local/include
INCL3 = ${HOME}/.local/include
LIB1 = ${HOME}/.local/lib
LIB2 = ${HOME}/.local/usr/lib
DEP  = petsc

install: Main.o Solver.o
	${LINKER} ${LFlags} Main.o Solver.o -l${DEP} -L${LIB1} -L${LIB2} -o${EXEC}

Main.o: Main.cpp
	${CC} ${CFlags} Main.cpp -I${INCL1} -I${INCL2} -I${INCL3}

Solver.o: Solver.cpp
	${CC} ${CFlags} Solver.cpp -I${INCL1} -I${INCL2} -I${INCL3}

clean: Main.o Solver.o
	${RM} Main.o Solver.o ${EXEC}
