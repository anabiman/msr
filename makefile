CC = mpiCC
RM = rm
LINKER = mpiCC
CFlags = -c -g -std=c++11 -fpic 
LFlags = -g
DEP  = petsc_real

SO = python/_msr.so
EXEC = msr.a

# Edit include + lib dirs below
INCL = -I/usr/lib/petscdir/3.10/include -I/usr/include/python3.7
LIB = /usr/lib/petscdir/3.10/lib

.PHONY: clean all

all: install python

install: main.o solver.o inout.o
	${LINKER} ${LFlags} main.o solver.o inout.o -l${DEP} -L${LIB} -o${EXEC}

python: solver.o inout.o interface_wrap.o
	${CC} -shared solver.o inout.o interface_wrap.o -l${DEP} -L${LIB} -o${SO}

src/interface_wrap.cpp:
	swig -python -c++ -o src/interface_wrap.cpp -outdir python python/interface.i

main.o: src/main.cpp
	${CC} ${CFlags} src/main.cpp ${INCL}

solver.o: src/solver.cpp
	${CC} ${CFlags} src/solver.cpp ${INCL}

inout.o: src/inout.cpp
	${CC} ${CFlags} src/inout.cpp ${INCL}

interface_wrap.o: src/interface_wrap.cpp
	${CC} ${CFlags} src/interface_wrap.cpp ${INCL}

clean:
	${RM} *.o
	${RM} ${EXEC}
	${RM} ${SO}
