# matrix multiplication with openmp intended as a Caliper test

CC=icpc
CC=g++
MPICC=mpicc
INC=
LIB=

all: mm mm_foo

mm: mm.c
	${CC} -g -o mm ${INC} mm.c -DORDER=3000 -fopenmp ${LIB}

mm_foo: mm_foo.c
	${CC} -O0 -g -o mm_foo ${INC} mm_foo.c -DORDER=3000 ${LIB} -fopenmp

mm_mpi: mm_mpi.c
	${MPICC} -O0 -g -o mm_mpi ${INC} mm_mpi.c -DORDER=1000 ${LIB} -fopenmp	

clean:
	rm -f mm mm_foo *.o
	# rm -rf MULTI__*

