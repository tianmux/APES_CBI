#include <parallel/algorithm>
#include <parallel/numeric>
const int OMP_NUM_THREADS = omp_get_max_threads()/2;

void update_M1(double* particles, double* M1, int nPar, int fill_count,int* Npar_survive,int* parFlag){
    double nParInver = 1.0/nPar;
    for(int i = 0;i<fill_count;++i){
        Npar_survive[i] = std::accumulate(parFlag+i*nPar,parFlag+(i+1)*nPar,0.0);
    }
    double* temppar = new double[nPar*fill_count*2];
    for (int i = 0;i<fill_count*nPar;++i){
        temppar[i] = particles[i]*parFlag[i];
        temppar[fill_count*nPar+i] = particles[fill_count*nPar+i]*parFlag[i];
    }
    for(int i = 0;i<fill_count;++i){
        M1[i] = std::accumulate(temppar+i*nPar,temppar+(i+1)*nPar,0.0)/Npar_survive[i];
        M1[fill_count+i] = std::accumulate(temppar+(fill_count+i)*nPar,temppar+(fill_count+i+1)*nPar,0.0)/Npar_survive[i];
    }
    delete[] temppar;
}

void update_M1_0(double* particles, double* M1, int nPar, int fill_count,int* Npar_survive,int* parFlag){
    double nParInver = 1.0/nPar;
    for(int i = 0;i<fill_count;++i){
        Npar_survive[i] = std::accumulate(parFlag+i*nPar,parFlag+(i+1)*nPar,0.0);
    }
    
    for(int i = 0;i<fill_count;++i){
        M1[i] = std::accumulate(particles+i*nPar,particles+(i+1)*nPar,0.0)*nParInver;
        M1[fill_count+i] = std::accumulate(particles+(fill_count+i)*nPar,particles+(fill_count+i+1)*nPar,0.0)*nParInver;
    }
}

void update_M2(double* particles, double* M1, double* M2, int nPar, int fill_count,int* Npar_survive,int* parFlag){
    // M2 store the sig_t_sq, sig_gamma_sq and sig_t_sig_gamma
    double nParInver = 1.0/nPar;
#pragma omp parallel for schedule(static,nPar) num_threads(OMP_NUM_THREADS)
    for(int i = 0;i<fill_count;++i){
        for(int j = 0;j<nPar;++j){
            M2[i] += (particles[i*nPar+j]-M1[i])*(particles[i*nPar+j]-M1[i])*parFlag[i*nPar+j];//time
            M2[fill_count+i] += (particles[fill_count*nPar+i*nPar+j]-M1[fill_count+i])*(particles[fill_count*nPar+i*nPar+j]-M1[fill_count+i])*parFlag[i*nPar+j]; // gamma
            M2[fill_count*2+i] += (particles[i*nPar+j]-M1[i])*(particles[fill_count*nPar+i*nPar+j]-M1[fill_count+i])*parFlag[i*nPar+j];//cross correlation
        }
        M2[i] = M2[i]/Npar_survive[i];
        M2[fill_count+i] = M2[fill_count+i]/Npar_survive[i];
        M2[fill_count*2+i] = M2[fill_count*2+i]/Npar_survive[i];
    }
    
}
