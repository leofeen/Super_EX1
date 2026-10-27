OMP_xl : 
	xlC_r -O3 -qstrict -qarch=pwr8 -qtune=pwr8 -qsmp=omp -I ./include -std=c++11 ./OMP.cpp ./source/solver.cpp -o OMP

OMP_gnu : 
	g++ -O3 -std=c++23 -fopenmp -I ./include -std=c++11 ./OMP.cpp ./source/solver.cpp -o OMP