for i in {100..2000..100}
do
	./seq_dgemm $i >> performance/seq_dgemm.txt
	./omp_dgemm $i >> performance/omp_dgemm.txt
done

