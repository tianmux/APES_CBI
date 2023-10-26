#ifndef INPUTPARA_H
#define INPUTPARA_H
#include <json/json.h>
#include <fstream>
#include <iostream>


// Class of input parameters
class inputPara{
public:
    double c = 299792458;
    double E0 = 938.2720813e6; //rest energy of proton in unit of eV
    double q = 1.6021766208e-19; // charge of proton
    int cT = 0; // current turn

    // general parameters
    std::map<std::string,unsigned int> generalPara;
   
    // ring parameters
    std::map<std::string,double> ringPara;
    // RF parameters
    std::map<std::string,std::vector<double>> rfPara;
    std::map<std::string,std::vector<double>> apPara;
    // Modulation
    std::map<std::string,std::vector<double>> modPara;

    // bunch parameters

    std::map<std::string,double> bunchPara;
    

    //beam pattern
    std::map<std::string,std::vector<int>> PatternPara;
    
    //ramp parameters
    std::map<std::string,std::vector<std::vector<double>>> RampPara;

    inputPara(){
        c = 299792458;
        E0 = 938.2720813e6; //rest energy of proton
        q = 1.6021766208e-19; // charge of proton
        cT = 0; // current turn
    
        // general parameters
        generalPara={
            {"type",1}, // 1 for electron, 0 for proton
            {"dynamicOn",1},
            {"n_dynamicOn",9000},
            {"n_Ek_damp_on",0},
            {"n_Ek_damp_off",3000},
            {"n_turns",1e8},
            {"turn_start",0},
            {"n_per_show",2},
            {"n_bunches",1},
            {"n_fill",10},// number of turn we start to fill the bunch into the bucket.
            {"n_q_ramp",2000},
            {"n_detune_start",40000}, // turn number when we start to ramp the frequency
            {"n_detune_ramp",80000}, // number of turns it takes to ramp the frequency.
            {"detune_slow_factor",10},
            {"n_detune_ramp_tot",100000}, // number of turns for second ramp.
            {"n_I_ramp_start",2000},// turn number when we start to ramp the Ig
            {"n_I_ramp_end",10000},// turn number when we end the ramp on Ig.
            {"mainRF",0}, // the main RF used for initial bunch distribution calculation.
            {"main_detune",0},
			{"Plot",0},
            {"init",0}, // whether this is a simulation from some previous result (0 means not from previous)
            {"write_break",1}, // whether we write the breaking point info to the disk.
            {"step_store",1},
            {"n_threads",18},
 			{"Slice",0},
            {"csv_out",0} // flag for the output type.
        };
       
        // ring parameters
        ringPara={
			{"R",610.1754}, // radius of the ring
			{"GMTSQ",552.25}, // gamma_t square
			{"Gamma",293},
			{"Ek",0},//ringPara["Gamma"]*E0},
			{"f0",0},//c/(2*M_PI*ringPara["R"])},
			{"T0",0},// 1/ringPara["f0"]},
			{"omg0",0},//2*M_PI*ringPara["f0"]},
			{"eta",0},//1/ringPara["GMTSQ"]-1/(ringPara["Gamma"]*ringPara["Gamma"])}
			{"nRF",1},
            {"nRF1",1},
			{"nRF2",1},
			{"nRFc",1},
            {"Prad",1e7}, // synchrotron radiation power.
            {"t_rad_long",1},
            {"Ek_damp",3.0},
            {"Ek_damp_always",5.0},
            {"nHOM",0}
        };

        // RF parameters
        rfPara={
            {"h",std::vector<double>(1)}, // the harmonic numbers of the RF, if there is only one then it essentially is the total number of buckets.
			{"RoQ",std::vector<double>(0)},
            {"Qs",std::vector<double>(0)},// This is the synchrotron tune, sqrt(h[0]*V[0]*eta*abs(cos(phis[0]))/(2*M_PI*Ek));
			{"QL",std::vector<double>(0)}, // this is the loaded Q of each mode. 
            {"Vref_I",std::vector<double>(0)},
			{"Vref_Q",std::vector<double>(0)} ,
            {"Iref_I",std::vector<double>(0)},
            {"Iref_Q", std::vector<double>(0)},
            {"I_I_ref_ini", std::vector<double>(0)},
            {"I_I_ref_final",std::vector<double>(0)},
            {"I_Q_ref_ini", std::vector<double>(0)},
            {"I_Q_ref_final",std::vector<double>(0)},
            {"n_ini_CBI",std::vector<double>(1)}, // number of initialized coupled bunch mode
            {"mu",std::vector<double>(1)}, // the mode index of initialized coupled bunch mode
            {"CBI_ini_amp",std::vector<double>(1)} // the amplitude of initial coupled bunch mode.
        };

        apPara={
            {"delay",std::vector<double>(1)}, // delay of the feedback.
            {"Update_Interval",std::vector<double>(1,1)},// the update interval in unit of time step, reflect the clock speed of the feedback control. 
			{"detune",std::vector<double>(1)},
            {"detune_ini",std::vector<double>(1)},
            {"detune_mid",std::vector<double>(1)}, // intermediate frequency where we switch tuning speed. 
            {"detune_final",std::vector<double>(1)},
            {"PA_cap",std::vector<double>(1)}, // the cap of the power amplifier output. in unit of Iref's
            {"nCav", std::vector<double>(1)},// number of cavities, for the radiation calculation. 
            {"epsilon_comb",std::vector<double>(0)},// the epsilon for the comb filter
            {"g_comb",std::vector<double>(0)}, // gain of the comb filter.
            {"n_fb_on",std::vector<double>(1)}
            
		};
        // Modulation
        modPara={
        {"Vepsilon",std::vector<double>(0)}, // modulation depth
        {"Mum",std::vector<double>(0)}, // modulation tune factor
        };
    
        // bunch parameters
    
        bunchPara={
            {"nBeam",1},
            {"beam_shift",0}, // this is the number of bucket between bunches in a train. 
            {"gap",0},
            {"Npar",1e6},
			{"NperBunch",2.7e9},
			{"N_bins",1e3},
            {"A",0.8},
            {"fill_step",1},
            {"siglong",1},
            {"epsln",bunchPara["A"]/6},//bunchPara["A"]/6},
            {"delta_hat",0},//sqrt(epsln)*sqrt(ringPara["omg0"]/(M_PI*ringPara["Ek"]))*std::pow((rfPara["h"][0]*rfPara["V"][0]*abs(cos(rfPara["phis"][0]))/(2*M_PI*ringPara["Ek"]*ringPara["eta"])),0.25)},
            {"t_hat",0},//bunchPara["delta_hat"]/rfPara["Qs"]*ringPara["eta"]/ringPara["omg0"]
        };

        // need to specify the beam patter, 
            // current idea is to specify number of train first, then specify the number of bunches and the gap between each train explicitly.
            // for example, 
            // "nTrain"=2;
            // "Pattern" = [1000,200,1000,200] 
            // means we have 2 trains, each has 1000 bunches in there and the gap between trains is 200 empty bucket.
        PatternPara={
            {"nTrain",std::vector<int>(1)},
            {"Pattern",std::vector<int>(2)},
        };
        RampPara = {
            {"nCurveType",std::vector<std::vector<double>>(1,std::vector<double>(1))},
            {"curveType",std::vector<std::vector<double>>(1,std::vector<double>(3))},// need at least three numbers to specify the ramp curve, 
                                                // first one represents the number of sections in this curve
                                                // the following pairs are the starting number of turn and the starting value of the ramped variable.
            {"idx_of_Type",std::vector<std::vector<double>>(1,std::vector<double>(1))}, //record the index of the bunch belong to specific ramping curve. 
        };
    };
	bool inRFPara(std::vector<std::string> V,std::vector<std::string> T, std::vector<std::string> TP,std::vector<std::string> Phi,std::string para){
		for (int i = 0;i<V.size();++i){
			if (V[i] == para|T[i] == para|TP[i] == para|Phi[i]==para){ 
				return true;
			}
		}
		return false;
	};
    int read(std::string path){
        // Object to write in file
        std::ifstream fin;
        fin.open(path, std::ios::in);
		std::vector<std::string> rfVLines;
		std::vector<std::string> rfTLines;
		std::vector<std::string> rfTPLines;
		std::vector<std::string> rfPhiLines;
        std::string temp;
        // Read the object's data in file
        while(std::getline(fin,temp)){
            if(temp.size()>1){
                std::stringstream ss(temp);
                std::istream_iterator<std::string> begin(ss);
                std::istream_iterator<std::string> end;
                std::vector<std::string> vstrings(begin, end);
                if(vstrings.size()>1){
                    if(vstrings[0]=="nRF"){
						int nRF = stoi(vstrings[1]);
						rfVLines.resize(nRF);
						rfTLines.resize(nRF);
						rfTPLines.resize(nRF);
						rfPhiLines.resize(nRF);

						for (int i = 0;i<nRF;++i){
							rfVLines[i] = "V"+std::to_string(i);
							rfTLines[i] = "TV"+std::to_string(i);
							rfTPLines[i] = "TP"+std::to_string(i);
							rfPhiLines[i] = "Phi"+std::to_string(i);
						}
                        rfPara["RoQ"].resize(nRF);
						rfPara["h"].resize(nRF);
						rfPara["QL"].resize(nRF);
                        rfPara["Vref_I"].resize(nRF);
                        rfPara["Vref_Q"].resize(nRF);
                        rfPara["Iref_I"].resize(nRF);
                        rfPara["Iref_Q"].resize(nRF);
                        rfPara["I_I_ref_ini"].resize(nRF);
                        rfPara["I_I_ref_final"].resize(nRF);
                        rfPara["I_Q_ref_ini"].resize(nRF);
                        rfPara["I_Q_ref_final"].resize(nRF);
                        apPara["fa"].resize(nRF);
                        apPara["phia"].resize(nRF);
                        apPara["n_ramp"].resize(nRF);
                        apPara["amplitude"].resize(nRF);
                        apPara["detune"].resize(nRF);
                        apPara["detune_ini"].resize(nRF);
                        apPara["detune_mid"].resize(nRF);
                        apPara["detune_final"].resize(nRF);
                        apPara["delay"].resize(nRF);
                        apPara["Update_Interval"].resize(nRF);
                        apPara["n_fb_on"].resize(nRF);
                        apPara["PA_cap"].resize(nRF);
                        apPara["nCav"].resize(nRF);
                        apPara["epsilon_comb"].resize(nRF);
                        apPara["g_comb"].resize(nRF);
                        apPara["gII"].resize(nRF);
                        apPara["gQQ"].resize(nRF);
                        apPara["gIQ"].resize(nRF);
                        apPara["gQI"].resize(nRF);
                        apPara["gIIi"].resize(nRF);
                        apPara["gQQi"].resize(nRF);
                        apPara["gIQi"].resize(nRF);
                        apPara["gQIi"].resize(nRF);
                    }
                    if(vstrings[0]=="nTrain"){
                        PatternPara["nTrain"][0]=stoi(vstrings[1]);
                        PatternPara["Pattern"].resize(stoi(vstrings[1])*3);
                    }

                    //std::cout<<vstrings[0]<<std::endl;
					if(inRFPara(rfVLines,rfTLines,rfTPLines,rfPhiLines,vstrings[0])){
						rfPara[vstrings[0]].resize(vstrings.size()-1);
						for(unsigned int i=0;i<vstrings.size()-1;++i){
                            rfPara[vstrings[0]][i]=stod(vstrings[i+1]); // store the ramping of the voltage.
                        }
					}
                    if(vstrings[0]=="n_ini_CBI"){
                        rfPara["n_ini_CBI"][0] = stoi(vstrings[1]);
                        rfPara["mu"].resize(rfPara["n_ini_CBI"][0]);
                        rfPara["CBI_ini_amp"].resize(rfPara["n_ini_CBI"][0]);
                    }

                    if(vstrings[0]=="mu"|vstrings[0]=="CBI_ini_amp"|vstrings[0]=="I_I_ref_ini"|vstrings[0]=="I_I_ref_final"|vstrings[0]=="I_Q_ref_ini"|vstrings[0]=="I_Q_ref_final"|vstrings[0]=="Iref_I"|vstrings[0]=="Iref_Q"|vstrings[0]=="RoQ"|vstrings[0]=="h"|vstrings[0]=="QL"|vstrings[0]=="Vref_I"|vstrings[0]=="Vref_Q"){
                        for(unsigned int i=0;i<rfPara[vstrings[0]].size();++i){
                            rfPara[vstrings[0]][i]=stod(vstrings[i+1]);
                        }
                    }

                    if(vstrings[0]=="csv_out"|vstrings[0]=="n_dynamicOn"|vstrings[0]=="dynamicOn"|vstrings[0]=="n_I_ramp_start"|vstrings[0]=="n_I_ramp_end"|vstrings[0]=="n_threads"|vstrings[0]=="main_detune"|vstrings[0]=="detune_slow_factor"|vstrings[0]=="n_detune_ramp_tot"|vstrings[0]=="mainRF"|vstrings[0]=="type"|vstrings[0]=="n_detune_start"|vstrings[0]=="n_detune_ramp"|vstrings[0]=="n_q_start"|vstrings[0]=="n_q_ramp"|vstrings[0]=="turn_start"|vstrings[0]=="init"|vstrings[0]=="step_store"|vstrings[0]=="Slice"|vstrings[0]=="Plot"|vstrings[0]=="n_turns"|vstrings[0]=="n_per_show"|vstrings[0]=="n_bunches"|vstrings[0]=="n_fill"|vstrings[0]=="turn_btw_fill"){
                        generalPara[vstrings[0]]=int(stod(vstrings[1]));
                    }

                    if(vstrings[0]=="Ek_damp_always"|vstrings[0]=="Ek_damp"|vstrings[0]=="nRF1"|vstrings[0]=="nRFc"|vstrings[0]=="nRF2"|vstrings[0]=="t_rad_long"|vstrings[0]=="Prad"|vstrings[0]=="nHOM"|vstrings[0]=="R"|vstrings[0]=="GMTSQ"|vstrings[0]=="Gamma"|vstrings[0]=="nRF"){
                        ringPara[vstrings[0]]=stod(vstrings[1]);
                    }

                    if(vstrings[0]=="Vepsilon"|vstrings[0]=="Mum"||vstrings[0]=="Tesp"||vstrings[0]=="Tmu"){
                        modPara[vstrings[0]].resize(vstrings.size()-1);
						for(unsigned int i=0;i<vstrings.size()-1;++i){
							modPara[vstrings[0]][i]=stod(vstrings[i+1]);
                        }
                    }

                    if(vstrings[0]=="beam_shift"|vstrings[0]=="nBeam"|vstrings[0]=="NperBunch"|vstrings[0]=="siglong"|vstrings[0]=="fill_step"|vstrings[0]=="Npar"|vstrings[0]=="A"|vstrings[0]=="N_bins"|vstrings[0]=="N_bnches_p_train"|vstrings[0]=="xmin"|vstrings[0]=="xmax"|vstrings[0]=="ymin"|vstrings[0]=="ymax"){
                        bunchPara[vstrings[0]]=stod(vstrings[1]);
                    }

					if(vstrings[0]=="Update_Interval"|vstrings[0]=="epsilon_comb"|vstrings[0]=="g_comb"|vstrings[0]=="nCav"|vstrings[0]=="PA_cap"|vstrings[0]=="detune_mid"|vstrings[0]=="detune_ini"|vstrings[0]=="detune_final"|vstrings[0]=="n_fb_on"|vstrings[0]=="gIIi"|vstrings[0]=="gQQi"|vstrings[0]=="gIQi"|vstrings[0]=="gQIi"|vstrings[0]=="gII"|vstrings[0]=="gQQ"|vstrings[0]=="gIQ"|vstrings[0]=="gQI"|vstrings[0]=="delay"|vstrings[0]=="fa"|vstrings[0]=="phia"|vstrings[0]=="amplitude"|vstrings[0]=="detune"|vstrings[0]=="n_ramp"){
                        for(unsigned int i=0;i<apPara[vstrings[0]].size();++i){
                            apPara[vstrings[0]][i]=stod(vstrings[i+1]);
                        }
                    }
                    if(vstrings[0]=="nTrain"){
                        for(unsigned int i=0;i<PatternPara[vstrings[0]].size();++i){                      
                            PatternPara[vstrings[0]][i]=stod(vstrings[i+1]);
                        }
                    }
                    if(vstrings[0]=="Pattern"){
                        for(unsigned int i=0;i<PatternPara[vstrings[0]].size();++i){
                            PatternPara[vstrings[0]][i]=stod(vstrings[i+1]);
                        }
                    }
                    if(vstrings[0]=="nCurveType"){
                        RampPara[vstrings[0]][0][0]=stod(vstrings[1]);
                        RampPara["curveType"].resize(RampPara[vstrings[0]][0][0]);
                        RampPara["idx_of_Type"].resize(RampPara[vstrings[0]][0][0]);
                    }
                    if(vstrings[0].substr(0,9)=="curveType"){
                        std::cout<<vstrings[0]<<std::endl;
                        int id = stod(vstrings[0].substr(9,vstrings[0].length()-9));
                        RampPara["curveType"][id].resize(stod(vstrings[1])*2+3);
                        for(unsigned int i=0;i<RampPara["curveType"][id].size();++i){
                            RampPara["curveType"][id][i]=stod(vstrings[i+1]);
                        }
                    }

                    if(vstrings[0].substr(0,11)=="idx_of_Type"){
                        std::cout<<vstrings[0]<<std::endl;
                        int id = stod(vstrings[0].substr(11,vstrings[0].length()-11));
                        std::cout<<"index list ID:  "<<id<<','<<std::endl;
                        if(id==0){// default curve type
                            std::cout<<"Setting up the default index list..."<<std::endl;
                            RampPara["idx_of_Type"][id].resize(generalPara["n_bunches"]);
                            std::cout<<"Total number of bunches : "<<generalPara["n_bunches"]<<std::endl;
                            for(int i = 0;i<generalPara["n_bunches"];++i){
                                RampPara["idx_of_Type"][id][i] = i;
                            }
                        }
                        else{
                            RampPara["idx_of_Type"][id].resize(vstrings.size()-1);
                            for(int i = 0;i<RampPara["idx_of_Type"][0].size();i++){
                                std::cout<<RampPara["idx_of_Type"][0][i]<<',';
                            }
                            std::cout<<std::endl;
                            for(unsigned int i=0;i<RampPara["idx_of_Type"][id].size();++i){
                                std::cout<<i<<','<<std::endl;
                                std::cout<<stod(vstrings[i+1])<<','<<std::endl;
                                RampPara["idx_of_Type"][id][i]=stod(vstrings[i+1]);
                                std::cout<<"Removing the id from default list."<<std::endl;
                                RampPara["idx_of_Type"][0].erase(std::find(RampPara["idx_of_Type"][0].begin(),RampPara["idx_of_Type"][0].end(),RampPara["idx_of_Type"][id][i]));
                            }
                            std::cout<<"Size of new default list : "<<RampPara["idx_of_Type"][0].size()<<std::endl;
                        }
                    }
                }
            }
        }
        ringPara["Ek"]=ringPara["Gamma"]*E0;
        std::cout<<"1."<<std::endl;
        ringPara["f0"]=(c*sqrt(1-1/ringPara["Gamma"]/ringPara["Gamma"]))/(2*M_PI*ringPara["R"]);
        std::cout<<"2."<<std::endl;
        ringPara["T0"]=1/ringPara["f0"];
        std::cout<<"3."<<std::endl;
        ringPara["omg0"]=2*M_PI*ringPara["f0"];
        std::cout<<"4"<<std::endl;
        ringPara["eta"]=1/ringPara["GMTSQ"]-1/(ringPara["Gamma"]*ringPara["Gamma"]);
        std::cout<<"5."<<std::endl;

        rfPara["Qs"].push_back(sqrt(rfPara["h"][0]*rfPara["Vref_I"][0]*ringPara["eta"]/(2*M_PI*ringPara["Ek"])));
        std::cout<<"6."<<std::endl;
		bunchPara["epsln"]=bunchPara["A"]/6.0;
        std::cout<<"7."<<std::endl;
        bunchPara["delta_hat"]=sqrt(bunchPara["epsln"])*sqrt(ringPara["omg0"]/(M_PI*ringPara["Ek"]))*std::pow((rfPara["h"][0]*rfPara["Vref_I"][0]/(2*M_PI*ringPara["Ek"]*ringPara["eta"])),0.25);
        std::cout<<"8."<<std::endl;
        bunchPara["t_hat"]=bunchPara["delta_hat"]/rfPara["Qs"][0]*ringPara["eta"]/ringPara["omg0"];
		return 0;
    };
    int printout(){
        std::cout<<"General Parameters:"<<std::endl;
        for(auto& x:generalPara){
            std::cout<<x.first<<"="<<x.second<<std::endl;
        }
        std::cout<<"Ring Parameters:"<<std::endl;
        for(auto& x:ringPara){
            std::cout<<x.first<<"="<<x.second<<std::endl;
        }
        std::cout<<"Rf Parameters:"<<std::endl;
        for(auto& x:rfPara){
            std::cout<<x.first;
            for(auto&data:x.second){
                std::cout<<":"<<data;
            }
            std::cout<<std::endl;
        }
		for(auto& x:apPara){
            std::cout<<x.first;
            for(auto&data:x.second){
                std::cout<<":"<<data;
            }
            std::cout<<std::endl;
        }
        std::cout<<"Modulation Parameters:"<<std::endl;
        for(auto& x:modPara){
            std::cout<<x.first;
            for(auto&data:x.second){
                std::cout<<":"<<data;
            }
            std::cout<<std::endl;
        }
        std::cout<<"Bunch Parameters:"<<std::endl;
        for(auto& x:bunchPara){
            std::cout<<x.first<<"="<<x.second<<std::endl;
        }
        std::cout<<"Pattern Parameters:"<<std::endl;
        for(auto& x:PatternPara){
            std::cout<<x.first;
            for(auto&data:x.second){
                std::cout<<":"<<data;
            }
            std::cout<<std::endl;
        }
        std::cout<<"Ramping Parameters:"<<std::endl;
        for(auto& x:RampPara){
            std::cout<<x.first;
            for(auto&data:x.second){
                std::cout<<std::endl;
                for(int i=0;i<data.size();++i)
                    std::cout<<":"<<data[i];
            }
            std::cout<<std::endl;
        }
        return 0;
    };
};

#endif