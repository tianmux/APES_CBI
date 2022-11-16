
#include <parallel/algorithm>
#include <parallel/numeric>

const int OMP_NUM_THREADS = omp_get_max_threads()/2;

void update_index_neighbor(double* bunch, double* weight, double dt, double Trev, int fill_count, 
                            int nPar, int totBin, int shift, int nBunch,int nBin,int shift_beam,int nBucket,
                            int* index, int* neighbor, int* parFlag){
    double dtInver = 1/dt;
    int totN = fill_count*nPar;
    // first totN data in bunch array is 'time' info of the particles.
#pragma omp parallel for schedule(static,nPar) num_threads(OMP_NUM_THREADS)
#pragma ivdep
    for(int i = 0;i<totN;++i){
        // first locate the particle on time grid. becareful that the particle time coordinates can be out of the bound
        //index[i] = (bunch[i] - Trev*fmod(bunch[i],Trev))*dtInver;
        index[i] = (totBin+int((fmod(bunch[i],Trev))*dtInver))%totBin;
#if 1
        if(index[i]<0|index[i]>totBin){
            //std::cout<<"Particle "<<i<<" is lost. "<<std::endl;
            parFlag[i] = 0;
        }
#endif
        // the we need to decide the weight of the particle on the grid

        weight[i] = fmod(bunch[i],Trev)-(dt*index[i]+dt*0.5);
        //std::cout<<weight[i]<<','<<dt<<','<<bunch[i]<<','<<Trev<<','<<fmod(bunch[i],Trev)<<','<<(dt*index[i]+dt*0.5)<<','<<index[i]<<std::endl;;
        int x = weight[i]>0?1:-1;
        neighbor[i] = (index[i] + x)%totBin;
        index[i] += shift;
        neighbor[i] += shift;
    }
}