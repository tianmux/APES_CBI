void getIQ(double* data, double* sintable, double* costable, double* tempI, double* tempQ, double* I, double* Q, 
            int delay, int nSmp, 
            int i, int nRF, int stride,int nBeam){
    int j = (i-delay-1)%nSmp;
    double nSmpInver = 1.0/nSmp;
#pragma vector always
#pragma ivdep
    for(int rf = 0;rf<nRF;++rf){
        int vidx = (i-delay)*stride+1+2*nBeam+3*rf;
        I[rf] = I[rf] - (tempI[j*nRF+rf]-data[vidx]*sintable[i*nRF+rf])*nSmpInver*2;
        Q[rf] = Q[rf] - (tempQ[j*nRF+rf]-data[vidx]*costable[i*nRF+rf])*nSmpInver*2;
        tempI[j*nRF+rf] = data[vidx]*sintable[i*nRF+rf];
        tempQ[j*nRF+rf] = data[vidx]*costable[i*nRF+rf];
    }
}