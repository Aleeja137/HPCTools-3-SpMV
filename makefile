CC=gcc
CFLAGS=-O2 -Wall -Wextra
LDLIBS=-lgslcblas

SRC=spmv.c my_dense.c my_sparse.c timer.c
OBJ=$(SRC:.c=.o)


spmv: $(OBJ)
	@echo $(OBJ)
	$(CC) $(CFLAGS) $(LDLIBS) -o $@ $^
clean:
	$(RM) $(OBJ) *~

cleanall:
	$(RM) $(OBJ) spmv *~