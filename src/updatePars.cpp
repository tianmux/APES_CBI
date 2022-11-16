#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>
#include "const.h"

void updateParticle(int type_of_particle, double* particle_coords, double* data, double* weight, double* M1,
                    double Gamma0, double eta, double dt, double Trev, double Ek_damp,double Urad0,
                    int h,int fill_step,int fill_count, int nPar, int nRF, int nHOM,int nBeam,int nRF1, int nRFc, int nRF2,int stride,int qovermp_on,double* nCav,
                    int* index, int* neighbor, int* parFlag){

    int nTot = fill_count*nPar;
    //int stride = 3+3*nRF+2*nHOM; // time, Ib[i-1], Ib[i], V0, Ig0, U0,
    double dtInver = 1/dt;
    double Gamma0Inver = 1/Gamma0;
    double nParInv = 1/nPar; 
    double qoverm;
    switch(type_of_particle)
    {
        case 0:
            qoverm = qovermp;
            break;
        case 2:
            qoverm = qovermAu;
            break;
        default:
            qoverm = qovermp;
    }
    // gamma of first beam
    int j = 0;
    for(int j = 0;j<nRF1+nRFc;++j){
#pragma omp parallel for //schedule(static,nPar) num_threads(OMP_NUM_THREADS)
#pragma ivdep
        for(int i = 0;i<nTot/nBeam;++i){
            int par_idx = index[i];
            int par_neighbor = neighbor[i];
            double V = 0;
            V = data[par_idx*stride+1+2*nBeam+j*3]+(data[par_neighbor*stride+1+2*nBeam+j*3]-data[par_idx*stride+1+2*nBeam+j*3])*dtInver*fabs(weight[i]);
            particle_coords[nTot+i] += qoverm*(V-Urad0*nCav[j])*c_inver*c_inver*qovermp_on;// Urad0 represent the unaccounted HOM power
        }
    }

    // gamma of second beam
    for(int j = nRF1;j<nRF1+nRFc+nRF2;++j){
#pragma omp parallel for //schedule(static,nPar) num_threads(OMP_NUM_THREADS)
#pragma ivdep
        for(int i = nTot/nBeam;i<nTot;++i){
            int par_idx = index[i];
            int par_neighbor = neighbor[i];
            double V = 0;
            V = data[par_idx*stride+1+2*nBeam+j*3]+(data[par_neighbor*stride+1+2*nBeam+j*3]-data[par_idx*stride+1+2*nBeam+j*3])*dtInver*fabs(weight[i]);
            particle_coords[nTot+i] = particle_coords[nTot+i]-qoverm*(V-Urad0*nCav[j])*c_inver*c_inver*qovermp_on; // '-' sign in front of qovermp means this is going in opposite direction.
        }
    }

    // artificial dampping
#pragma omp parallel for //schedule(static,nPar) num_threads(OMP_NUM_THREADS)
#pragma ivdep 
    for(int i = 0;i<nTot;++i){
        particle_coords[nTot+i] -= (M1[fill_count+int(i*nParInv)]-Gamma0)*Ek_damp*qovermp_on; // damping term
    }

    // t
#pragma omp parallel for //schedule(static,nPar) num_threads(OMP_NUM_THREADS) 
    for(int i = 0;i<nTot;++i){
        int iBunch = i/nPar;
        double tstart = Trev/h*iBunch*fill_step;
        double tend = Trev/h*(iBunch*fill_step+1);
        particle_coords[i] += Trev/(1/(eta*(particle_coords[nTot+i]*Gamma0Inver-1))-1)*qovermp_on;
        if(particle_coords[i]<tstart | particle_coords[i]>tend){
            //std::cout<<particle_coords[i]<<','<<tstart<<','<<tend<<std::endl;
            //std::cout<<i<<','<<iBunch<<std::endl;
            //exit(-1);
            parFlag[i] = 0;
        }
        //particle_coords[i] = fmod(particle_coords[i]+Trev,Trev);
    }
    
    //update_M1(particle_coords, M1, nPar, fill_count);
}
void updateParticle_e(double* particle_coords, double* data, double* weight, double* M1,double* rand_table,
                    double Gamma0, double eta, double dt, double Trev, double Ek_damp,double Urad0,double coeff_decay_long, double coeff_excite_long,
                    int fill_count, int nPar, int nRF, int nHOM,int nBeam,int nRF1, int nRFc, int nRF2,int stride,int qoverme_on,double* nCav,
                    int* index, int* neighbor){

    int nTot = fill_count*nPar;
    //int stride = 3+3*nRF+2*nHOM; // time, Ib[i-1], Ib[i], V0, Ig0, U0,
    double dtInver = 1/dt;
    double Gamma0Inver = 1/Gamma0;
    double nParInv = 1/nPar; 
    
    // gamma of first beam
    for(int j = 0;j<nRF1+nRFc;++j){
#pragma omp parallel for //schedule(static,nPar) num_threads(OMP_NUM_THREADS)
#pragma ivdep
        for(int i = 0;i<nTot/nBeam;++i){
            int par_idx = index[i];
            int par_neighbor = neighbor[i];
            //std::cout<<par_idx<<','<<par_neighbor<<std::endl;
            double V = 0;
            V = data[par_idx*stride+1+2*nBeam+j*3]+(data[par_neighbor*stride+1+2*nBeam+j*3]-data[par_idx*stride+1+2*nBeam+j*3])*dtInver*fabs(weight[i]);
            particle_coords[nTot+i] += qoverme*(V-Urad0*nCav[j])*c_inver*c_inver*qoverme_on;//radiation
        }
    }

    // artificial dampping
#pragma omp parallel for //schedule(static,nPar) num_threads(OMP_NUM_THREADS)
#pragma ivdep 
    for(int i = 0;i<nTot;++i){
        particle_coords[nTot+i] -= (M1[fill_count+int(i*nParInv)]-Gamma0)*Ek_damp*qoverme_on; // damping term
       
    }

// minic the sychrotron radiation dampping and quantum excitation 
#pragma omp parallel for //schedule(static,nPar) num_threads(OMP_NUM_THREADS)
#pragma ivdep 
    for(int i = 0;i<nTot;++i){
        particle_coords[nTot+i] -= qoverme_on*(coeff_decay_long*(particle_coords[nTot+i]-Gamma0)+coeff_excite_long*(2*double(rand_table[i])-1));
    }

    // tQs0
#pragma omp parallel for //schedule(static,nPar) num_threads(OMP_NUM_THREADS) 
    for(int i = 0;i<nTot;++i){
        particle_coords[i] += Trev/(1/(eta*(particle_coords[nTot+i]*Gamma0Inver-1))-1)*qoverme_on;
        particle_coords[i] = fmod(particle_coords[i]+Trev,Trev);
    }
}

void rampCharge(int nCurveType, std::vector<std::vector<double>> ramp_type,std::vector<std::vector<int>> idx_of_Type,std::vector<double> &qpb,double qpb0,int turn,int nBunches){
    // look at the bunch index,
    // see what ramp curve this bunch should follow,
    // then take the turn number, 
    // then figure out what value of qpb should this bunch have now.
#pragma omp parallel for //schedule(static,nPar) num_threads(OMP_NUM_THREADS)
#pragma ivdep 
    for(int i = 0;i<nBunches;++i){
        for(int j = 0;j<nCurveType;++j){
            auto pos=std::find(idx_of_Type[j].begin(),idx_of_Type[j].end(),i);
            if(pos!=idx_of_Type[j].end()){// if the bunch belongs to this type of ramp
                //std::cout<<"Bunch "<<i<<" follows type "<<j<<" curve."<<std::endl;
                int nSec = ramp_type[j][0]; // number of section in this type of ramp.
                for(int k = 0;k<nSec;++k){
                    //for(int m = 0;m<ramp_type[j].size();m++){
                    //   std::cout<<ramp_type[j][m]<<',';
                    //}
                    //std::cout<<std::endl;
                    //std::cout<<ramp_type[j][1+k*2]<<',';
                    //std::cout<<ramp_type[j][1+k*2+2]<<',';
                    //std::cout<<turn<<std::endl;;
                    if(ramp_type[j][1+k*2]<turn && ramp_type[j][1+(k+1)*2]>=turn){
                        //std::cout<<"Gradient is "<<(ramp_type[j][1+(k+1)*2+1]-ramp_type[j][1+k*2+1])/(ramp_type[j][1+(k+1)*2]-ramp_type[j][1+k*2])<<std::endl;
                        //std::cout<<qpb0*(ramp_type[j][1+(k+1)*2+1]-ramp_type[j][1+k*2+1])/(ramp_type[j][1+(k+1)*2]-ramp_type[j][1+k*2])*(turn-ramp_type[j][1+k*2])<<std::endl;
                        qpb[i]=qpb0*(ramp_type[j][1+k*2+1]+(ramp_type[j][1+(k+1)*2+1]-ramp_type[j][1+k*2+1])/(ramp_type[j][1+(k+1)*2]-ramp_type[j][1+k*2])*(turn-ramp_type[j][1+k*2]));
                        //std::cout<<i<<": "<<qpb[i]<<std::endl;
                        break;
                    }
                }
                break;
            }

        }

    }
    

}