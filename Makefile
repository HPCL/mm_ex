# matrix multiplication with openmp intended as a Caliper test

CC=icpc
CC=g++
MPICC=mpicxx
INC=
LIB=

USE_CALI=True

ifdef ORDER
	INC+=-DORDER=${ORDER}
endif


ifdef USE_CALI 
  INC += -I${CALIPER_DIR}/include
  INC += -DUSE_CALI
  LIB  += -L${CALIPER_DIR}/lib64 
  LIB  += -lcaliper -lcaliper-mpi
endif


all: mm mm_foo

mm: mm.c
	${CC} -g -o mm ${INC} mm.c -DORDER=3000 -fopenmp ${LIB}

mm_foo: mm_foo.c
	${CC} -O0 -g -o mm_foo ${INC} mm_foo.c -DORDER=3000 ${LIB} -fopenmp

mm_mpi: mm_mpi.c
	${MPICC} -O2 -g -o mm_mpi ${INC} mm_mpi.c ${LIB} -fopenmp	

clean:
	rm -f mm mm_foo mm_mpi *.o
	# rm -rf MULTI__*

