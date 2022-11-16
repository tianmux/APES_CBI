#include<vector>
#include<complex>
#include<iostream>
#include<omp.h>
#include<cmath>
void update_Vb(std::vector<std::complex<double>> &Vb, std::vector<double> t, std::vector<double> Vadd, double wrf, double alpha){
    std::complex<double> V0 = Vb[Vb.size()-1];
    double t0 = t[t.size()-1];
    std::complex<double> j(0,1);
    std::complex<double> iwrf = j*wrf-alpha;
    Vb[0] = V0*exp((j*wrf-alpha)*(t[0]))+Vadd[0];
    for(int i=1;i<Vb.size();i++){
        Vb[i] = Vb[i-1]*std::exp((iwrf)*(t[i]))+Vadd[i];
    }
}

int main(int argc, char* argv[]){
    auto t_start = omp_get_wtime(); 
    std::vector<std::complex<double>> Vb;
    std::vector<double> t;
    std::vector<double> Vadd;

    const double wrf = 2*M_PI*100;
    const double alpha = 1e-10;

    for( int count = 0; count < argc; count++ )
         std::cout << "  argv[" << count << "]   "
                << argv[count] << "\n";
    const int N = atoi(argv[1]);
    for(int i=0;i<N;i++){
        Vb.push_back(0);
        t.push_back(i);
        Vadd.push_back(0);
    }
    Vb[N-1] = std::complex<double>(1,0);
    
    update_Vb(Vb,t,Vadd,wrf,alpha);
    auto t_end = omp_get_wtime(); 
    std::cout<<"Total running time : "<<(t_end-t_start) <<" [s]. "<<std::endl;
    std::cout<<"Average time for one turn : "<< (t_end-t_start)*1000000/N<<" [us]. "<<std::endl;
    std::cout<<Vb[0]<<std::endl;
    std::cout<<Vb[N-1]<<std::endl;

    
    return 0;
}