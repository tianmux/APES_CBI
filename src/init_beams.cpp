#include <random>
#include <vector>
#include <iostream>
#include "const.h"

void init_bunch(double* bunch,double t0,double gamma0,double sigt,double siggamma,int Npar, int nBunch,int bunch_idx,int nBeam){
    /* So far only the Gaussian distribution is supported */
    std::default_random_engine generator(1); // Normal distribution generator
    generator.seed(0);
    std::normal_distribution<double> dist0(t0, sigt); // normal distribution in t
    std::normal_distribution<double> dist1(0, siggamma); // normal distribution in gamma    
    
    for(int i = 0;i<Npar;++i){
        bunch[bunch_idx*Npar+i] = dist0(generator);
        
        while(abs(bunch[bunch_idx*Npar+i]-t0)>3*sigt){
            //std::cout<<"re_throw"<<std::endl;
            bunch[bunch_idx*Npar+i] = dist0(generator);
        }
    }
    for(int i = 0 ;i<Npar;++i){
        bunch[nBunch*nBeam*Npar+bunch_idx*Npar+i] = (dist1(generator)+1)*gamma0;//0.0*cos(2*pi*bunch_idx/120)+gamma0;
        
        while(abs(bunch[nBunch*nBeam*Npar+bunch_idx*Npar+i]-gamma0)>3*siggamma*gamma0){
            bunch[nBunch*nBeam*Npar+bunch_idx*Npar+i] = (dist1(generator)+1)*gamma0;
        }
        
    }

    bunch[bunch_idx*Npar] = t0; // set the first particle to the reference particle.
    bunch[nBunch*nBeam*Npar+bunch_idx*Npar] = gamma0;
    //bunch[nBunch*Npar+bunch_idx] = gamma0+1;
}

void init_tempIQ(double* data, double* sintable, double* costable, double* tempI, double* tempQ, 
                double* I, double* Q, 
                int* delay, int* nSmp, int i,int nRF,int stride){
// initialize the temp arrays that hold the intermediate results for I Q calculation.    
    for(int j = 0;j<nRF;++j){    
        int idx_temp = 0;
        for(int idx = i-delay[0]-nSmp[0]+1;idx<i-delay[0]+1;++idx){
        tempI[idx_temp*nRF+j] = data[idx*stride+3+j*3]*sintable[idx*nRF+j];
        tempQ[idx_temp*nRF+j] = data[idx*stride+3+j*3]*costable[idx*nRF+j];
        idx_temp++;
        }
        for(int idx = 0;idx<nSmp[0];++idx){
            I[j] += tempI[idx*nRF+j];
            Q[j] += tempQ[idx*nRF+j];
        }
        I[j] = I[j]/nSmp[0]*2;
        Q[j] = Q[j]/nSmp[0]*2;
    }
}

void get_init_t0s_backup(double* t0s,double dT, double t0, double Trev, int nBeam, int nBunch, int shift,int nBucket,std::vector<int> &mode, std::vector<double> &amp,int n_ini_CBI = 0){
// get the initial time coordinates of each bunches stored in array "t0s", dT is the bunch distance, t0 is the start t
// nBeam is the number of beams, can be one or two
// nBunch is the number of bunches in each beam
// shift is the relative shift between two bunches 
// 'mode' is the mu number for coupled bunch mode initially excited
// 'amp' is the amplitude of the coupled bunch mode, can excite only one mode now.
    std::cout<<"shift is : "<<shift<<std::endl;
    for(int i = 0;i<nBeam;++i){
        for(int j = 0;j<nBunch;++j){
            t0s[i*nBunch+j] = fmod(t0+(i*shift)%nBucket*t0*2+j*dT+i*t0/2,Trev);
            // apply initial coupled bunch mode if any
            for(int k = 0;k<n_ini_CBI;++k){
                t0s[i*nBunch+j] += amp[k]*sin(mode[k]*2*pi*1/Trev*t0s[i*nBunch+j]);
            }
            
        }

        std::cout<<"Time of first bunch: "<<t0s[i*nBunch+0]<<std::endl;
    }
}
void get_init_t0s(double* t0s,double dT, double t0, double Trev, 
                int nBeam, int nBunch, int shift,
                int nBucket,std::vector<int> &mode, std::vector<double> &amp,int n_ini_CBI,
                int nTrain,std::vector<int> &pattern){
// get the initial time coordinates of each bunches stored in array "t0s", dT is the bunch distance, t0 is the start t
// nBeam is the number of beams, can be one or two
// nBunch is the number of bunches in each beam
// shift is the relative shift between two bunches 

// 'mode' is the mu number for coupled bunch mode initially excited
// 'amp' is the amplitude of the coupled bunch mode, can excite only one mode now.

// nTrain is the number of bunch trains.
// pattern stores the information of bunch pattern, in the form of nBunch1, fill_step1, nGap1, nBunch2, fill_step2, nGap2, ..., in unit of RF bucket. 
    for(int i = 0; i<pattern.size();++i){
        std::cout<<pattern[i]<<","<<std::endl;
    }
    for(int i = 0;i<nBeam;++i){
        int bunch_index = 0;
        int bucket_index = 0;
        for(int j = 0;j<nTrain;++j){
            for(int n = 0 ; n<pattern[j*3];++n){
                std::cout<<"bunchId: "<<bunch_index<<" ; "<<"bucketID: "<<bucket_index<<std::endl;
                //t0s[i*nBunch+j] = fmod(t0+(i*shift)%nBucket*t0*2+j*dT+i*t0/2,Trev);
                //this is for one beam case only.
                t0s[i*nBunch+bunch_index] = t0+bucket_index*dT;
                // apply initial coupled bunch mode if any, might have some problem with non-uniform fill.
                for(int k = 0;k<n_ini_CBI;++k){
                    t0s[i*nBunch+bunch_index] += amp[k]*sin(mode[k]*2*pi*1/Trev*t0s[i*nBunch+bunch_index]);
                }
                bunch_index += 1;
                bucket_index += pattern[j*3+1];
            }
            bucket_index += pattern[j*3+2];
        }
        std::cout<<"Sanity check: nBunch consistent?"<<nBunch-bunch_index<<std::endl;
        std::cout<<"Time of first bunch: "<<t0s[i*nBunch+0]<<std::endl;
    }
}
void get_init_gamma0s(double* gamma0s,double* t0s,double dT, double t0, double Trev, double Gamma0, int nBeam, int nBunch, int shift,int nBucket,std::vector<int> &mode, std::vector<double> &amp,int n_ini_CBI = 0){
// get the initial gamma coordinates of each bunches stored in array "gamma0s", dT is the bunch distance, t0 is the start t
// nBeam is the number of beams, can be one or two
// nBunch is the number of bunches in each beam
// shift is the relative shift between two bunches 
    std::cout<<"shift is : "<<shift<<std::endl;
    for(int i = 0;i<nBeam;++i){
        for(int j = 0;j<nBunch;++j){
            gamma0s[i*nBunch+j] = Gamma0;
            for(int k = 0;k<n_ini_CBI;++k){
                gamma0s[i*nBunch+j] += amp[k]*cos(mode[k]*2*pi/Trev*t0s[i*nBunch+j]);
            }
            
        }
    }
}


void current_rampping(int turn, int n_I_ramp_start, int n_I_ramp_end, int nRF, double dt, 
                        double* I_I_ref, double* I_Q_ref, double* I_I_ref_ini, double* I_I_ref_final, double* I_Q_ref_ini, double* I_Q_ref_final){
    double rampInv = 1.0/double(n_I_ramp_end-n_I_ramp_start);
    #pragma ivdep 
    for(int i = 0;i<nRF;++i){
        I_I_ref[i] = (turn-n_I_ramp_start)*(I_I_ref_final[i]-I_I_ref_ini[i])*rampInv+I_I_ref_ini[i];
        I_Q_ref[i] = (turn-n_I_ramp_start)*(I_Q_ref_final[i]-I_Q_ref_ini[i])*rampInv+I_Q_ref_ini[i];
    }
}

void detune_rampping(int turn, int n_detune_start, int n_detune_ramp,int nRF, double dt,
            double* omegac, double* omegarf, double* QL, double* C, double* R, double* Leff,double* k,
            double* detune, double* detune_ini, double* detune_final){
    double rampInv = 1.0/double(n_detune_ramp-n_detune_start);
#pragma ivdep 
    for(int i = 0;i<nRF;++i){
        detune[i] = (turn-n_detune_start)*(detune_final[i]-detune_ini[i])*rampInv+detune_ini[i];
        omegac[i] = omegarf[i]+detune[i]*2*pi;
        C[i] = QL[i]/R[i]/omegac[i];//[2] = {QL[0]/R[0]/omegac[0],QL[1]/R[1]/omegac[1]};//9.2755373928124199e-11;
        //Leff[i] = dt*dt/(1-cos(dt*omegac[i]))/2/C[i];//1/(C[i]*omegac[i]*omegac[i]);//
        //Leff[i] = 1/(omegac[i]*omegac[i]*C[i]);
        Leff[i] = dt*dt/(2-2*cos(dt*omegac[i])-dt/R[i])/C[i];
        k[i*5] = 1-dt/R[i]/C[i];
        k[i*5+1] = dt/C[i];
        k[i*5+2] = -dt/C[i]/2;
        k[i*5+3] = k[i*5+2];
        k[i*5+4] = -dt/C[i]/Leff[i];
    }
}
