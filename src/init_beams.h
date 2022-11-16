#include <vector>

void init_bunch(double* bunch,double t0,double gamma0,double sigt,double siggamma,int Npar, int nBunch,int bunch_idx,int nBeam);
void init_tempIQ(double* data, double* sintable, double* costable, double* tempI, double* tempQ, 
                double* I, double* Q, 
                int* delay, int* nSmp, int i,int nRF,int stride);
void get_init_t0s_backup(double* t0s,double dT, double t0, double Trev, int nBeam, int nBunch, int shift,int nBucket,std::vector<int> &mode, std::vector<double> &amp,int n_ini_CBI = 0);
void get_init_t0s(double* t0s,double dT, double t0, double Trev, 
                int nBeam, int nBunch, int shift,
                int nBucket,std::vector<int> &mode, std::vector<double> &amp,int n_ini_CBI,
                int nTrain,std::vector<int> &pattern);
void get_init_gamma0s(double* gamma0s,double* t0s,double dT, double t0, double Trev, double Gamma0, int nBeam, int nBunch, int shift,int nBucket,std::vector<int> &mode, std::vector<double> &amp,int n_ini_CBI = 0);
void current_rampping(int turn, int n_I_ramp_start, int n_I_ramp_end, int nRF, double dt, 
                        double* I_I_ref, double* I_Q_ref, double* I_I_ref_ini, double* I_I_ref_final, double* I_Q_ref_ini, double* I_Q_ref_final);

void detune_rampping(int turn, int n_detune_start, int n_detune_ramp,int nRF, double dt,
            double* omegac, double* omegarf, double* QL, double* C, double* R, double* Leff,double* k,
            double* detune, double* detune_ini, double* detune_final);