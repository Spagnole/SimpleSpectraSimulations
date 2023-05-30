#include "ReadInputData.cc"

//this function fills the simulated peak intensities
//this function requires the output from the ReadDecayScheme(), GetEfficiency() and GetPeakWidth()
//scaling paramter is used to fix for using relative efficiency and/or relative intensities
//for absolute efficiency and absolute intensities just use the default value
void FillSimulation(double scale_int = 1.0){

	if( AssignedTransition.size() == 0 ){
		cout << "No data in decay list!\nHave you read in list of gamma-rays?" << endl;
		return;
	}

	const int Nbins = hBkgr->GetXaxis()->GetNbins();
	double x_low = hBkgr->GetXaxis()->GetBinLowEdge(1);
	double x_max = hBkgr->GetXaxis()->GetBinUpEdge(Nbins);
	hSimPeaks = new TH1D("hSimPeaks","Simulated Peaks",Nbins, x_low, x_max);
	double peak_int;
	double energy;
	double width;
	for(int i = 0; i < AssignedTransition.size(); i++){
		energy = get<1>(AssignedTransition.at(i)); 										//get gamma-ray energy
		peak_int = get<3>(AssignedTransition.at(i)) * gEff->Eval(energy) * scale_int; 	//calculate peak intensity
		width = fWidth->Eval(energy); 													//get peak width
		for(int j = 0; j < peak_int; j++){
			hSimPeaks->Fill( gRandom->Gaus(energy,width) ); 							//fill spectrum using simple gaussian distribution
		}
	}
	hSimPeaks->SetLineColor(kRed);
}

//this function fills the simulated single-escape peak intensities
//this function requires the output from the ReadEscapePeaks() function
//!!!USE SAME SCALING THAT WAS USED BY FillSimulation()
void FillEscapePeaks(double scale_int = 1.0){

	const int Nbins = hBkgr->GetXaxis()->GetNbins();
	double x_low = hBkgr->GetXaxis()->GetBinLowEdge(1);
	double x_max = hBkgr->GetXaxis()->GetBinUpEdge(Nbins);
	hEscPeaks = new TH1D("hEscPeaks","Simulated Escape Peaks",Nbins, x_low, x_max);
	hEscPeaks->SetLineColor(kGreen+2);

	double peak_int;
	double energy;
	double width;
	for(int i = 0; i < AssignedTransition.size(); i++){
		if( get<1>(AssignedTransition.at(i)) < 1500. ) continue;
		energy = get<1>(AssignedTransition.at(i)); 															//get gamma-ray energy
		peak_int = get<3>(AssignedTransition.at(i))*gEff->Eval(energy)*fEscPeak->Eval(energy) * scale_int;	//calculate escape-peak intensity
		width = 1.1*fWidth->Eval(energy);																	//get peak width, S.E.P width is increased by 10 %
		for(int j = 0; j < peak_int; j++){
			hEscPeaks->Fill( gRandom->Gaus(energy-511.,width) );											//fill spectrum using simple gaussian distribution
		}
	}
}


//this function is used to produce the full simulated spectrum
//this function reads in the simulated Full energy peaks, single escape peaks and the background extracted from the real spectrum
//this function creates three new spectra
//One adds the full-energy peaks, single esc, peaks and the background, this is the full simulated spectrum
//Another adds the full energy peaks to the background
//The third adds the escape peaks to the background
void BuildSimuledSpectra(TH1D *hSimPeaks, TH1D *hBG, TH1D *hEscPeak){


	const int Nbins = hBkgr->GetXaxis()->GetNbins();
	double x_low = hBkgr->GetXaxis()->GetBinLowEdge(1);
	double x_max = hBkgr->GetXaxis()->GetBinUpEdge(Nbins);
	hFullSim = new TH1D("hFullSim","Full Sim Spectrum",Nbins, x_low, x_max);
	hSimPeaks_Bkgr = new TH1D("hSimPeaks_Bkgr","Sim Peaks on Bkgr",Nbins, x_low, x_max);
	hEscPeaks_Bkgr = new TH1D("hEscPeaks_Bkgr","Sim Esc, Peaks on Bkgr",Nbins, x_low, x_max);


	hFullSim->SetLineColor(kRed);
	hFullSim->Add(hBG);
	hFullSim->Add(hSimPeaks);
	hFullSim->Add(hEscPeak);

	hSimPeaks_Bkgr->SetLineColor(kGreen+2);
	hSimPeaks_Bkgr->Add(hBG);
	hSimPeaks_Bkgr->Add(hSimPeaks);
		
	hEscPeaks_Bkgr->SetLineColor(6);
	hEscPeaks_Bkgr->Add(hBG);
	hEscPeaks_Bkgr->Add(hEscPeak);

	hRealSpectra->SetFillColor(kBlue);
	hRealSpectra->SetFillStyle(3003);
		
	hRealSpectra->Draw("hist");
	hSimPeaks_Bkgr->Draw("histsame");
	hEscPeaks_Bkgr->Draw("histsame");		
	hFullSim->Draw("histsame");
	hBG->Draw("histsame");
}



