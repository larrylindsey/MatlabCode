CC=g++
CFLAGS= -Wall -O3
INCLUDE= -I${EM_CODE_DIR}/lib/libnpy/include
LIBNPY= ${EM_CODE_DIR}/lib/libnpy/lib/libnpy.a
BINDIR= ${EMROOT}/bin
OBJECTS=testnpy.o

all: testnpy

testnpy: testnpy.o ${LIBNPY}
	 ${CC} ${CFLAGS} -o ${BINDIR}/$@ $< ${LIBNPY} 

# rule for object-source dependency
%.o: %.cpp
	 ${CC} ${CFLAGS} ${INCLUDE} -c $< -o $@
