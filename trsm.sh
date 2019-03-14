for i in {100..2000..100}
do
	./seq_dtrsm $i U >> performance/seq_dtrsm.txt
	./omp_dtrsm $i U >> performance/omp_dtrsm.txt
done

