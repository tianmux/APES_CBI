tmake: 
	dpcpp -qopenmp ./src/main.cpp -o main.o -qopenmp -std=c++20 -O3 -g -xCORE-AVX2 -I /opt/intel/advisor/include -qmkl=parallel -c
	dpcpp -qopenmp ./src/getIQ.cpp -o getIQ.o -qopenmp -std=c++20 -O3 -g -xCORE-AVX2 -I /opt/intel/advisor/include -qmkl=parallel -c
	dpcpp -qopenmp ./src/updateVI.cpp -o updateVI.o -qopenmp -std=c++20 -O3 -g -xCORE-AVX2 -I /opt/intel/advisor/include -qmkl=parallel -c
	dpcpp -qopenmp ./src/update_index_neighbor.cpp -o update_index_neighbor.o -qopenmp -std=c++20 -O3 -g -xCORE-AVX2 -I /opt/intel/advisor/include -qmkl=parallel -c
	dpcpp -qopenmp ./src/updateMs.cpp -o updateMs.o -qopenmp -std=c++20 -O3 -g -xCORE-AVX2 -I /opt/intel/advisor/include -qmkl=parallel -c
	dpcpp -qopenmp ./src/updatePars.cpp -o updatePars.o -qopenmp -std=c++20 -O3 -g -xCORE-AVX2 -I /opt/intel/advisor/include -qmkl=parallel -c
	dpcpp -qopenmp ./src/init_beams.cpp -o init_beams.o -qopenmp -std=c++20 -O3 -g -xCORE-AVX2 -I /opt/intel/advisor/include -qmkl=parallel -c
	dpcpp -o APESICPX main.o getIQ.o updateVI.o update_index_neighbor.o updateMs.o updatePars.o init_beams.o -qopenmp -O3 -g -xCORE-AVX2 	
