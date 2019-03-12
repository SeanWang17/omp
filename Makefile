CC = icc
CCFLAGS = -mkl -O2 -openmp
ALL= omp_dgemm seq_dgemm omp_dtrsm seq_dtrsm
all: $(ALL)
.PHONY : all

omp_dgemm: omp_dgemm.c
	$(CC) $< -o $@ $(CCFLAGS) 

seq_dgemm: seq_dgemm.c
	$(CC) $< -o $@ $(CCFLAGS) 

omp_dtrsm: omp_dtrsm.c
	$(CC) $< -o $@ $(CCFLAGS) 

seq_dtrsm: seq_dtrsm.c
	$(CC) $< -o $@ $(CCFLAGS) 

clean:
	rm -f $(ALL) 
