# matrix multiplication with openmp intended as a Caliper test

CC=icpc
INC=-I${CALIPER_DIR}/include
LIB=-lcaliper
ORDER=1000

all: original

original: mm_original.c
	${CC} -O0 -g -o original ${INC} mm_original.c -DORDER=${ORDER} ${LIB} -fopenmp

swapped_loops: mm_swapped_loops.c
	${CC} -O0 -g -o swapped_loops ${INC} mm_swapped_loops.c -DORDER=${ORDER} ${LIB} -fopenmp

transpose: mm_transpose.c
	${CC} -O0 -g -o transpose ${INC} mm_transpose.c -DORDER=${ORDER} ${LIB} -fopenmp

callpaths: mm_callpaths.c
	${CC} -O0 -g -o callpaths ${INC} mm_callpaths.c -DORDER=${ORDER} ${LIB} -fopenmp

clean:
	rm -f mm mm_foo *.o *.cali *.json
	# rm -rf MULTI__*

