#include <vector>
#include <cmath>
#include "getIQ.h"
#include <iostream>

void updateVI(double *data,
                double* k, 
                int nData_tot,
                double dt, 
                double* I_I_ref,
                double* I_Q_ref,
                double* V_I_ref,
                double* V_Q_ref,
                double* PA_cap,
                double* R,
                double* gii,
                double* gqq,
                double* gii_i,
                double* gqq_i,
                double* I,
                double* Q,
                double* sintable,
                double* costable,
                double* tempI,
                double* tempQ,
                std::vector<std::vector<double>>& Sp,
                std::vector<std::vector<double>>& S,
                std::vector<double> epsilon_comb,
                std::vector<double> g_comb,
                std::vector<double>& errorI_p,
                std::vector<double>& errorQ_p,
                std::vector<double>& errorI_i,
                std::vector<double>& errorQ_i,
                int* delay, 
                int* Update_Interval,
                int* nSmp,
                int nRF, int nRF1,int nRFc, int nRF2,int stride,int nBeam){ // k is the constant parameter vector
    double* Rinv = new double[nRF];
    double* I_buffer = new double[nRF];
    double* Q_buffer = new double[nRF]; // used to store the I and Q "Update_Interval" time steps before 
    int iprv = 0;
    for(int i = 0;i<nRF;++i){
        Rinv[i] = 1/R[i];
    }

    // Update the V,I of every time step.
    for(int i = delay[0]+nSmp[0];i<nData_tot;++i){
        //first calculate the previous index
        iprv = i-1;
        // then update the voltage
        // the meanning of the terms are V(i-1), Ig(i-1), Ib(i), Ib(i-1),U(i-1)
#pragma ivdep       
        for(int j = 0;j<nRF1;++j){
            //data[i*stride+1+2*nBeam+j*3] = data[iprv*stride+1+2*nBeam+j*3]*k[j*5]+(data[iprv*stride+1+2*nBeam+j*3+1]+data[i*stride+1+2*nBeam+j*3+1])/2*k[j*5+1]+data[iprv*stride+1]*k[j*5+2]+data[iprv*stride+2]*k[j*5+3]+data[iprv*stride+1+2*nBeam+j*3+2]*k[j*5+4];
            //data[i*stride+1+2*nBeam+j*3] = data[iprv*stride+1+2*nBeam+j*3]*k[j*5]+(data[iprv*stride+1+2*nBeam+j*3+1]*0.5+data[i*stride+1+2*nBeam+j*3+1]*0.5)*k[j*5+1]+data[i*stride+1]*k[j*5+2]+data[i*stride+2]*k[j*5+3]+data[iprv*stride+1+2*nBeam+j*3+2]*k[j*5+4];
            data[i*stride+1+2*nBeam+j*3] = data[iprv*stride+1+2*nBeam+j*3]*k[j*5]+data[iprv*stride+1+2*nBeam+j*3+1]*k[j*5+1]+data[i*stride+1]*k[j*5+2]+data[i*stride+2]*k[j*5+3]+data[iprv*stride+1+2*nBeam+j*3+2]*k[j*5+4];
            //data[i*stride+1+2*nBeam+j*3] = data[iprv*stride+1+2*nBeam+j*3]*k[j*5]+data[i*stride+1+2*nBeam+j*3+1]*k[j*5+1]+data[i*stride+1]*k[j*5+2]+data[i*stride+2]*k[j*5+3]+data[iprv*stride+1+2*nBeam+j*3+2]*k[j*5+4];
            data[i*stride+1+2*nBeam+j*3+2] = data[iprv*stride+1+2*nBeam+j*3+2]+dt*data[i*stride+1+2*nBeam+j*3];
        }
#pragma ivdep
        for(int j = nRF1;j<nRF1+nRFc;++j){
            data[i*stride+1+2*nBeam+j*3] = data[iprv*stride+1+2*nBeam+j*3]*k[j*5]+data[iprv*stride+1+2*nBeam+j*3+1]*k[j*5+1]+(data[i*stride+1]+data[i*stride+3])*k[j*5+2]+(data[i*stride+2]+data[i*stride+4])*k[j*5+3]+data[iprv*stride+1+2*nBeam+j*3+2]*k[j*5+4];
            data[i*stride+1+2*nBeam+j*3+2] = data[iprv*stride+1+2*nBeam+j*3+2]+dt*data[i*stride+1+2*nBeam+j*3];
        }
#pragma ivdep
        for(int j = nRF1+nRFc;j<nRF1+nRFc+nRF2;++j){
            data[i*stride+1+2*nBeam+j*3] = data[iprv*stride+1+2*nBeam+j*3]*k[j*5]+data[iprv*stride+1+2*nBeam+j*3+1]*k[j*5+1]+data[i*stride+3]*k[j*5+2]+data[i*stride+4]*k[j*5+3]+data[iprv*stride+1+2*nBeam+j*3+2]*k[j*5+4];
            data[i*stride+1+2*nBeam+j*3+2] = data[iprv*stride+1+2*nBeam+j*3+2]+dt*data[i*stride+1+2*nBeam+j*3];
        }
        // then we need to get the I, Q to update current
        getIQ(data,sintable, costable, tempI, tempQ, I,Q,delay[0],nSmp[0],i,nRF,stride,nBeam);

        for(int j = 0; j<nRF; ++j){
            S[j][i-nSmp[0]] = (1-epsilon_comb[j])*Sp[j][i-nSmp[0]]+epsilon_comb[j]*data[i*stride+1+2*nBeam+j*3];
            Sp[j][i-nSmp[0]] = S[j][i-nSmp[0]];
        }
        // update Ig for all frequency. 
#pragma ivdep
        for(int j = 0;j<nRF;++j){
            // direct IQ feedback
            // check if we need to apply the feedback based on the Update_Interval parameter
            if (i%Update_Interval[j]==0){
                I_buffer[j] = I[j];
                Q_buffer[j] = Q[j];
            }

            // proportional 
            errorI_p[j] = (V_I_ref[j]-I_buffer[j]);
            errorQ_p[j] = (V_Q_ref[j]-Q_buffer[j]);
                                

            // integral part
            errorI_i[j] += errorI_p[j]*dt*Update_Interval[j];
            errorQ_i[j] += errorQ_p[j]*dt*Update_Interval[j];
            double Ig_I_correction = errorI_p[j]*Rinv[j]*gii[j] + errorI_i[j]*Rinv[j]*gii_i[j];
            double Ig_Q_correction = errorQ_p[j]*Rinv[j]*gqq[j] + errorQ_i[j]*Rinv[j]*gqq_i[j];
            data[i*stride+1+2*nBeam+j*3+1] = (I_I_ref[j]+Ig_I_correction)*sintable[i*nRF+j]+(I_Q_ref[j]+Ig_Q_correction)*costable[i*nRF+j];
            // comb

            //data[i*stride+1+2*nBeam+j*3+1] += g_comb[j]*Rinv[j]*(V_I_ref[j]*sintable[i*nRF+j]+V_Q_ref[j]*costable[i*nRF+j]-S[j][i-nSmp[0]-delay[0]]);
            data[i*stride+1+2*nBeam+j*3+1] += g_comb[j]*(I_I_ref[j]*sintable[i*nRF+j]+I_Q_ref[j]*costable[i*nRF+j]-Rinv[j]*(S[j][i-nSmp[0]-delay[0]]));
            //data[i*stride+1+2*nBeam+j*3+1] += g_comb[j]*Rinv[j]*S[j][i-nSmp[0]-delay[0]];

        }
    }
    delete[] Rinv;    
}


void updateIb(double* particle_coords, double* data, double* weight, 
            double dt, double* q, std::vector<double> &qpb,
            int fill_count, int nPar, int totBin, int shift, int stride,int nBeam,
            int* index, int* neighbor,int* parFlag){
    double dtInver = 1/dt;
    int nBunchPerBeam = int(fill_count/nBeam);
    // first reset the Ibs

    for(int j = 0;j<nBeam;++j){
#pragma vector always
#pragma ivdep
        for(int i = 0; i<totBin; ++i){
            data[i*stride+j*2+1] = 0;
            data[i*stride+j*2+2] = 0;
        }
    }

    for(int j = 0;j<nBeam;++j){
        int beam_start_idx = j*nBunchPerBeam*nPar;
        for(int i = beam_start_idx; i<beam_start_idx+1*nBunchPerBeam*nPar ; ++i){
            if(i/nPar==0){
                //std::cout<<index[i]<<std::endl;
            }
            int bunch_idx = int(i/nPar);
            //std::cout<<"index of bunch: "<<bunch_idx<<std::endl;
            //std::cout<<"qpb : "<<qpb[beam_start_idx+bunch_idx]<<std::endl;
            //std::cout<<"q : "<<q[j]<<std::endl;
            
            data[(index[i])*stride+j*2+1] += qpb[beam_start_idx+bunch_idx]*(1-fabs(weight[i]*dtInver))*dtInver*parFlag[i];
            data[(neighbor[i])*stride+j*2+1] += qpb[beam_start_idx+bunch_idx]*fabs(weight[i]*dtInver)*dtInver*parFlag[i];
            data[(index[i]+1)*stride+j*2+2] += qpb[beam_start_idx+bunch_idx]*(1-fabs(weight[i]*dtInver))*dtInver*parFlag[i];
            data[(neighbor[i]+1)*stride+j*2+2] += qpb[beam_start_idx+bunch_idx]*fabs(weight[i]*dtInver)*dtInver*parFlag[i];
        }
    }
}
void updateIb0(double* particle_coords, double* data, double* weight, 
            double dt, double* q, 
            int fill_count, int nPar, int totBin, int shift, int stride,int nBeam,
            int* index, int* neighbor,int* parFlag){
    double dtInver = 1/dt;
    int nBunchPerBeam = int(fill_count/nBeam);
    // first reset the Ibs

    for(int j = 0;j<nBeam;++j){
#pragma vector always
#pragma ivdep
        for(int i = 0; i<totBin; ++i){
            data[i*stride+j*2+1] = 0;
            data[i*stride+j*2+2] = 0;
        }
    }

    for(int j = 0;j<nBeam;++j){
        //std::cout<<q[j]<<std::endl;
        for(int i = j*nBunchPerBeam*nPar; i<(j+1)*nBunchPerBeam*nPar ; ++i){
            if(i/nPar==0){
                //std::cout<<index[i]<<std::endl;
            }
            data[(index[i])*stride+j*2+1] += q[j]*(1-fabs(weight[i]*dtInver))*dtInver*parFlag[i];
            data[(neighbor[i])*stride+j*2+1] += q[j]*fabs(weight[i]*dtInver)*dtInver*parFlag[i];
            data[(index[i]+1)*stride+j*2+2] += q[j]*(1-fabs(weight[i]*dtInver))*dtInver*parFlag[i];
            data[(neighbor[i]+1)*stride+j*2+2] += q[j]*fabs(weight[i]*dtInver)*dtInver*parFlag[i];
        }
    }
}
