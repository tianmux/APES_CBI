tmake: 
	dpcpp -qopenmp ./update_Vb.cpp -o main.o -qopenmp -std=c++20 -O3 -g -xCORE-AVX2 -I /opt/intel/advisor/include -qmkl=parallel -c
	
	dpcpp -o Cython_tests main.o -qopenmp -O3 -g -xCORE-AVX2 	
