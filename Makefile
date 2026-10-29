OMP_xl : 
	xlC_r -O3 -qstrict -qarch=pwr8 -qtune=pwr8 -qsmp=omp -I ./include -std=c++11 ./OMP.cpp ./source/solver.cpp -o OMP

xl : 
	xlC_r -O3 -qstrict -qarch=pwr8 -qtune=pwr8 -I ./include -std=c++11 ./OMP.cpp ./source/solver.cpp -o SEQ