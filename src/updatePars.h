void updateParticle(int type_of_particle, double* particle_coords, double* data, double* weight, double* M1,
                    double Gamma0, double eta, double dt, double Trev, double Ek_damp,double Urad0,
                    int h,int fill_step,int fill_count, int nPar, int nRF, int nHOM,int nBeam,int nRF1, int nRFc, int nRF2,int stride,int qovermp_on,double* nCav,
                    int* index, int* neighbor, int* parFlag);
void updateParticle_e(double* particle_coords, double* data, double* weight, double* M1,double* rand_table,
                    double Gamma0, double eta, double dt, double Trev, double Ek_damp,double Urad0,double coeff_decay_long, double coeff_excite_long,
                    int fill_count, int nPar, int nRF, int nHOM,int nBeam,int nRF1, int nRFc, int nRF2,int stride,int qovermp_on,double* nCav,
                    int* index, int* neighbor);
void rampCharge(int nCurveType, std::vector<std::vector<double>> ramp_type,std::vector<std::vector<int>> idx_of_Type,std::vector<double> &qpb,double qpb0,int turn,int nBunches);
