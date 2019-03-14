for i in {5..20}
do
	./omp_dgemm 1000 $i >> performance/gemm_weak.txt
	./omp_dtrsm 1000 1000 U $i >> performance/trsm_weak.txt
done
