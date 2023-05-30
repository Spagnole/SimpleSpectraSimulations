#include "ReadInputData.cc"

bool PrintReducedCoincList = false;

//this vector is a list of levels, this vector is populated by the GetLevelList() function.
//The components of this vector are: Level Index ---> Level Energy ---> Level Population
vector<tuple<int, double, double>> LevelEnergy;

//list of gamma-gamma coincidence intensities
//The components of this vector are: Gamma 1 ---> Gamma 2 ---> Intensity of coincidence
vector<tuple<double,double,double>> gg_coinc; 

TH2D *sim_gg_mat;
TH1D *SimMatProj;
TH1D *hEscPeak;
TH1D *hFullSimSpectrum;


//this function gets a list of levels from the 'AssignedTransition' vector from the ReadInputData.cc file
//this function requires that ReadDecayScheme() function has been utilized
void GetLevelList(){

	if( AssignedTransition.size() == 0 ){
		cout << "No data in decay list!\nHave you read in list of gamma-rays?" << endl;
		return;
	}

	double prev_energy = -1;
	int counter = 0;
	for(int i = 0; i < AssignedTransition.size(); i++){
		if( get<0>(AssignedTransition.at(i)) != prev_energy ){
			LevelEnergy.push_back( make_tuple(counter, get<0>(AssignedTransition.at(i)), 0) );
			counter++;
		}
		prev_energy = get<0>(AssignedTransition.at(i));
	}
	
	for(int i = 0; i < LevelEnergy.size(); i++){
		if(PrintLevelList) cout << get<0>(LevelEnergy.at(i)) << "\t" << get<1>(LevelEnergy.at(i)) <<endl;
	}
}

//this function is used to correct the 'AssignedTransition' list
//this correction is required if your final level energies do not exactly match the level energies in the 'LevelEnergy' vector

void FixTransitionList(){

	if( AssignedTransition.size() == 0 ){
		cout << "No data in decay list!\nHave you read in list of gamma-rays?" << endl;
		return;
	}

	double min_diff;
	double difference;
	double new_level_energy;
	for(int i = 0; i <  AssignedTransition.size(); i++){
		min_diff = 10000;
		for(int j = 0; j < LevelEnergy.size(); j++){
			difference = TMath::Abs( get<2>(AssignedTransition.at(i)) - get<1>(LevelEnergy.at(j)) );
			if( difference < min_diff ){
				new_level_energy = get<1>(LevelEnergy.at(j));
				min_diff = difference;
			}
		}
		get<2>(AssignedTransition.at(i)) = new_level_energy;		
	}
}


//this function calculates the level population of each state and the gamma-ray branching ratio of each transition
//this is function is required to properly calculate the expected gamma-gamma coincidence intensity
void CalcLevelFeedingAndGammaBR(){
	
	double sum_temp;
	for(int i = 0; i < LevelEnergy.size(); i++){
		sum_temp = 0;
		for(int j = 0; j <  AssignedTransition.size(); j++){
			if(  get<1>(LevelEnergy.at(i)) == get<0>(AssignedTransition.at(j))){
				sum_temp += get<3>(AssignedTransition.at(j));
			}
		}
		get<2>(LevelEnergy.at(i)) = sum_temp;
	}
	
	//calculating gamma-ray branching ratios
	for(int i = 0; i <  AssignedTransition.size(); i++){
		if( i == 0 ) continue;
		for(int j = 0; j < LevelEnergy.size(); j++){
			if(  get<1>(LevelEnergy.at(j)) == get<0>(AssignedTransition.at(i))){
				get<4>(AssignedTransition.at(i)) = get<2>(LevelEnergy.at(j));
				get<5>(AssignedTransition.at(i)) = get<3>(AssignedTransition.at(i)) / get<2>(LevelEnergy.at(j));
			}
		}
	}
	
	//calculate level population
	for(int i = 0; i <  AssignedTransition.size(); i++){
		if( i == 0 ) continue;		
		for(int j = 0; j < AssignedTransition.size(); j++){
			if(  get<2>(AssignedTransition.at(j)) == get<0>(AssignedTransition.at(i))){
				get<4>(AssignedTransition.at(i)) = get<4>(AssignedTransition.at(i)) - get<3>(AssignedTransition.at(j));
			}
		}
	}
	
	//calculating gamma-ray branching ratios
	for(int i = AssignedTransition.size()-1; i >= 0; i--){
		if(PrintCalcBR) cout << get<0>(AssignedTransition.at(i)) << "\t" << get<1>(AssignedTransition.at(i)) << "\t" << get<2>(AssignedTransition.at(i)) << "\t";
		if(PrintCalcBR) cout << get<3>(AssignedTransition.at(i)) << "\t" << get<4>(AssignedTransition.at(i)) << "\t" << get<5>(AssignedTransition.at(i)) << endl;
	}
}

//this function calculates the expected number of gamma-gamma coincidences
//this function uses a large number of nested loops (Should be replaced by some recursive loop?)
//if a gamma-gamma coincidence is separated by about 5 intermediate gamma-rays it will not be added to the list of coincidences
//if you require such coincidences ---> Add more nested loops
void FindCoincidences(){
		
	if( AssignedTransition.size() == 0 ){
		cout << "No data in decay list!\nHave you read in list of gamma-rays?" << endl;
		return;
	}

	double NDecays, Ngg[10];
	
	for(int i = AssignedTransition.size()-1; i >= 0 ; i--){
		if( get<2>(AssignedTransition.at(i))  == 0 ) continue;
		NDecays = get<4>(AssignedTransition.at(i))*get<5>(AssignedTransition.at(i));
		for(int j = i-1; j >= 0 ; j--){
			if( get<2>(AssignedTransition.at(i)) == get<0>(AssignedTransition.at(j))){
				Ngg[0] = NDecays*get<5>(AssignedTransition.at(j));
				gg_coinc.push_back(make_tuple(get<1>(AssignedTransition.at(i)), get<1>(AssignedTransition.at(j)), Ngg[0]));
				if(PrintFindingCoincidences) cout << get<0>(AssignedTransition.at(i)) << "\t" << get<1>(AssignedTransition.at(i)) << "\t" << get<2>(AssignedTransition.at(i)) << "\t";
				if(PrintFindingCoincidences) cout << get<1>(AssignedTransition.at(j)) << "\t" << get<2>(AssignedTransition.at(j)) << "\t" << Ngg[0] << endl;
				//cout << get<0>(AssignedTransition.at(i)) << "\t" << get<1>(AssignedTransition.at(i)) << "\t" << get<2>(AssignedTransition.at(i)) << endl;
				//cout << "\t" << get<0>(AssignedTransition.at(j)) << "\t" << get<1>(AssignedTransition.at(j)) << "\t" << get<2>(AssignedTransition.at(j)) << endl;
				if( get<2>(AssignedTransition.at(j))  == 0 ) continue;
				for(int k = j-1; k >= 0 ; k--){
					if( get<2>(AssignedTransition.at(j)) == get<0>(AssignedTransition.at(k))){
						Ngg[1] = Ngg[0]*get<5>(AssignedTransition.at(k));
						gg_coinc.push_back(make_tuple(get<1>(AssignedTransition.at(i)), get<1>(AssignedTransition.at(k)), Ngg[1]));
						gg_coinc.push_back(make_tuple(get<1>(AssignedTransition.at(j)), get<1>(AssignedTransition.at(k)), Ngg[1]));
						if(PrintFindingCoincidences) cout << "\t\t" << get<1>(AssignedTransition.at(k)) << "\t" << get<2>(AssignedTransition.at(k)) << "\t" << Ngg[1] << endl;
						//cout << "\t\t" << get<0>(AssignedTransition.at(k)) << "\t" << get<1>(AssignedTransition.at(k)) << "\t" << get<2>(AssignedTransition.at(k)) << endl;
						if( get<2>(AssignedTransition.at(k))  == 0 ) continue;
						for(int l = k-1; l >= 0 ; l--){
							if( get<2>(AssignedTransition.at(k)) == get<0>(AssignedTransition.at(l))){
								Ngg[2] = Ngg[1]*get<5>(AssignedTransition.at(l));
								gg_coinc.push_back(make_tuple(get<1>(AssignedTransition.at(i)), get<1>(AssignedTransition.at(l)), Ngg[2]));
								gg_coinc.push_back(make_tuple(get<1>(AssignedTransition.at(j)), get<1>(AssignedTransition.at(l)), Ngg[2]));
								gg_coinc.push_back(make_tuple(get<1>(AssignedTransition.at(k)), get<1>(AssignedTransition.at(l)), Ngg[2]));
								if(PrintFindingCoincidences) cout << "\t\t\t" << get<1>(AssignedTransition.at(l)) << "\t" << get<2>(AssignedTransition.at(l)) << "\t" << Ngg[2] << endl;
								//cout << "\t\t\t" << get<0>(AssignedTransition.at(l)) << "\t" << get<1>(AssignedTransition.at(l)) << "\t" << get<2>(AssignedTransition.at(l)) << endl;
								if( get<2>(AssignedTransition.at(l))  == 0 ) continue;
								for(int m = l-1; m >= 0 ; m--){
									if( get<2>(AssignedTransition.at(l)) == get<0>(AssignedTransition.at(m))){
										Ngg[3] = Ngg[2]*get<5>(AssignedTransition.at(m));
										gg_coinc.push_back(make_tuple(get<1>(AssignedTransition.at(i)), get<1>(AssignedTransition.at(m)), Ngg[3]));
										gg_coinc.push_back(make_tuple(get<1>(AssignedTransition.at(j)), get<1>(AssignedTransition.at(m)), Ngg[3]));
										gg_coinc.push_back(make_tuple(get<1>(AssignedTransition.at(k)), get<1>(AssignedTransition.at(m)), Ngg[3]));
										gg_coinc.push_back(make_tuple(get<1>(AssignedTransition.at(l)), get<1>(AssignedTransition.at(m)), Ngg[3]));
										if(PrintFindingCoincidences) cout << "\t\t\t\t" << get<0>(AssignedTransition.at(m)) << "\t" << get<1>(AssignedTransition.at(m)) << "\t" <<  Ngg[3] << endl;
										if( get<2>(AssignedTransition.at(m))  == 0 ) continue;
										for(int n = m-1; n >= 0 ; n--){
											if( get<2>(AssignedTransition.at(m)) == get<0>(AssignedTransition.at(n))){
												Ngg[4] = Ngg[3]*get<5>(AssignedTransition.at(n));
												gg_coinc.push_back(make_tuple(get<1>(AssignedTransition.at(i)), get<1>(AssignedTransition.at(n)), Ngg[4]));
												gg_coinc.push_back(make_tuple(get<1>(AssignedTransition.at(j)), get<1>(AssignedTransition.at(n)), Ngg[4]));
												gg_coinc.push_back(make_tuple(get<1>(AssignedTransition.at(k)), get<1>(AssignedTransition.at(n)), Ngg[4]));
												gg_coinc.push_back(make_tuple(get<1>(AssignedTransition.at(l)), get<1>(AssignedTransition.at(n)), Ngg[4]));
												gg_coinc.push_back(make_tuple(get<1>(AssignedTransition.at(m)), get<1>(AssignedTransition.at(n)), Ngg[4]));
												if(PrintFindingCoincidences) cout << "\t\t\t\t\t" << get<0>(AssignedTransition.at(n)) << "\t" << get<1>(AssignedTransition.at(n)) << "\t" << Ngg[4] << endl;
												if( get<2>(AssignedTransition.at(n))  == 0 ) continue;
												for(int p = n-1; p >= 0 ; p--){
													if( get<2>(AssignedTransition.at(n)) == get<0>(AssignedTransition.at(p))){
														Ngg[5] = Ngg[4]*get<5>(AssignedTransition.at(p));
														gg_coinc.push_back(make_tuple(get<1>(AssignedTransition.at(i)), get<1>(AssignedTransition.at(p)), Ngg[5]));
														gg_coinc.push_back(make_tuple(get<1>(AssignedTransition.at(j)), get<1>(AssignedTransition.at(p)), Ngg[5]));
														gg_coinc.push_back(make_tuple(get<1>(AssignedTransition.at(k)), get<1>(AssignedTransition.at(p)), Ngg[5]));
														gg_coinc.push_back(make_tuple(get<1>(AssignedTransition.at(l)), get<1>(AssignedTransition.at(p)), Ngg[5]));
														gg_coinc.push_back(make_tuple(get<1>(AssignedTransition.at(m)), get<1>(AssignedTransition.at(p)), Ngg[5]));
														gg_coinc.push_back(make_tuple(get<1>(AssignedTransition.at(n)), get<1>(AssignedTransition.at(p)), Ngg[5]));
														if( get<2>(AssignedTransition.at(p)) !=0){
															//cout  << "WARNING"<<endl;
															cout << "\t\t\t\t\t\t" << get<0>(AssignedTransition.at(n)) << "\t" << get<1>(AssignedTransition.at(n)) << "\t" <<Ngg[5] << endl;
														}
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}


//this function is used to reduce the size of the gamma-gamma coincidence list
//if a given coincidence is given a number of times, this function will reduce the entries of that coincidence to one instance but adds all of the intensities
//strictly speaking this function is not absolutely required for the tool to function! BUT YOU SHOULD STILL USE IT!
void ReduceGammaGammaList(){
	
	for(int i = 0; i < gg_coinc.size();i++){
		for(int j = 0; j < gg_coinc.size(); j++){
			if(i==j) continue;
			if( get<0>(gg_coinc.at(i)) == get<0>(gg_coinc.at(j)) && get<1>(gg_coinc.at(i)) == get<1>(gg_coinc.at(j)) ){
				if( PrintReducedCoincList ) cout << get<0>(gg_coinc.at(i)) << "\t" << get<1>(gg_coinc.at(i)) << "\t"  << get<2>(gg_coinc.at(i)) <<endl;
				if( PrintReducedCoincList ) cout <<"\t" << get<0>(gg_coinc.at(j)) << "\t" << get<1>(gg_coinc.at(j)) << "\t"  << get<2>(gg_coinc.at(j)) <<endl;
				get<2>(gg_coinc.at(i)) =  get<2>(gg_coinc.at(i)) + get<2>(gg_coinc.at(j));
				if( PrintReducedCoincList ) cout << "\t\t" << get<0>(gg_coinc.at(i)) << "\t" << get<1>(gg_coinc.at(i)) << "\t"  << get<2>(gg_coinc.at(i)) <<endl;
				gg_coinc.erase(gg_coinc.begin()+j);
				j=j-1;
			}
		}
	}
}

//this function fills a simulated gamma-gamma coincidence matrix
//user can provide specified binning of the matrix
//this function requires that the list of coincidences has been produced
//as well as reading in the gamma-ray efficiency and peak widths
//
void FillCoincMatrix(int NBins = 8100, double low = -0.5, double upp = 8099.5){
	
	sim_gg_mat = new TH2D("sim_gg_mat","Simulated g-g matrix",NBins,low,upp,NBins,low,upp);
	double gamma1,gamma2,sigma1,sigma2;
	double counts;
	for(int i = 0; i < gg_coinc.size(); i++){
		gamma1 = get<0>(gg_coinc.at(i));
		gamma2 = get<1>(gg_coinc.at(i));
		counts = (14./15) * get<2>(gg_coinc.at(i)) * gEff->Eval(gamma1) * gEff->Eval(gamma2); //user may require different normalisation if using relative efficiency curve
		sigma1 = fWidth->Eval(gamma1);
		sigma2 = fWidth->Eval(gamma2);
		for(int j = 0; j < counts; j++){
			sim_gg_mat->Fill( gRandom->Gaus(gamma1,sigma1), gRandom->Gaus(gamma2,sigma2));
			sim_gg_mat->Fill( gRandom->Gaus(gamma2,sigma2), gRandom->Gaus(gamma1,sigma1));
		}
	}
	
}


//this function gates on the simulated gamma-gamma matrix to produce a coincidence spectrum
//used should use the same gate they used for their experimental efficiency curve
void GateOnSimMat(int low, int up){
	
	SimMatProj  = sim_gg_mat->ProjectionY(Form("SimProj_%d_%d",low,up),low,up);
	SimMatProj->Draw("hist");
	SimMatProj->SetLineColor(kRed);
	SimMatProj->SetFillColor(kRed);
	SimMatProj->SetFillStyle(3000);
	
}

//this function fills the simulated escape peaks 
//this is done by using the output TH1D from the GateOnSimMat() function
void FillEscPeakSpec(TH1D *hSim){

	int Nbins = hSim->GetXaxis()->GetNbins();
	double low = hSim->GetXaxis()->GetBinLowEdge(1);
	double upp = hSim->GetXaxis()->GetBinUpEdge(Nbins);
	
	hEscPeak = new TH1D("hEscPeak","Simulated Esc. Peaks",Nbins, low, upp);
	hEscPeak->SetLineColor(6);
	hEscPeak->SetFillColor(6);
	hEscPeak->SetFillStyle(3000);
	double gamma1,sigma1;
	double counts, scale; 
	
	for(int i = 1; i <= Nbins; i++){
		gamma1 = hSim->GetBinCenter(i);
		if(gamma1 < 1500. ) continue;
		counts = hSim->GetBinContent(i);
		for(int j = 0; j < counts; j++){
			hEscPeak->Fill( gamma1-511. , fEscPeak->Eval(gamma1) );
		}
	}
}

//this function is identical to the function used by SimulateSingles.cc (Perhaps it should be put in the ReadInputData.cc file?)
//this function is used to produce the full simulated spectrum
//this function reads in the simulated Full energy peaks, single escape peaks and the background extracted from the real spectrum
//this function creates three new spectra
//One adds the full-energy peaks, single esc, peaks and the background, this is the full simulated spectrum
//Another adds the full energy peaks to the background
//The third adds the escape peaks to the background
void BuildSimuledSpectra(TH1D *hSimPeaks, TH1D *hBkgr, TH1D *hEscPeaks){

	const int Nbins = hBkgr->GetXaxis()->GetNbins();
	double x_low = hBkgr->GetXaxis()->GetBinLowEdge(1);
	double x_max = hBkgr->GetXaxis()->GetBinUpEdge(Nbins);
	hFullSim = new TH1D("hFullSim","Full Sim Spectrum",Nbins, x_low, x_max);
	hSimPeaks_Bkgr = new TH1D("hSimPeaks_Bkgr","Sim Peaks on Bkgr",Nbins, x_low, x_max);
	hEscPeaks_Bkgr = new TH1D("hEscPeaks_Bkgr","Sim Esc, Peaks on Bkgr",Nbins, x_low, x_max);


	hFullSim->SetLineColor(kRed);
	hFullSim->Add(hBkgr);
	hFullSim->Add(hSimPeaks);
	hFullSim->Add(hEscPeaks);

	hSimPeaks_Bkgr->SetLineColor(kGreen+2);
	hSimPeaks_Bkgr->Add(hBkgr);
	hSimPeaks_Bkgr->Add(hSimPeaks);
		
	hEscPeaks_Bkgr->SetLineColor(6);
	hEscPeaks_Bkgr->Add(hBkgr);
	hEscPeaks_Bkgr->Add(hEscPeaks);

	hRealSpectra->SetFillColor(kBlue);
	hRealSpectra->SetFillStyle(3003);
		
	hRealSpectra->Draw("hist");
	hSimPeaks_Bkgr->Draw("histsame");
	hEscPeaks_Bkgr->Draw("histsame");		
	hFullSim->Draw("histsame");
	hBkgr->Draw("histsame");
}
