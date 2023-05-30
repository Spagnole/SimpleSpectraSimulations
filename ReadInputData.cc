//these bool statements is simply used for debugging
//
bool PrintReadData = false;
bool PrintLevelList = false;
bool PrintCalcBR = false;
bool PrintFindingCoincidences = false;
bool PrintPeakWidths = false;
bool PrintEscPeakData = false;


//////////////////////////////////////////////////////////////
///// PLEASE NOT ALL ENERGIES SHOULD BE GIVEN IN keV!!! //////
//////////////////////////////////////////////////////////////


//This is for the list of transitions from your decay scheme
//there are six components to this vector
//Initial level energy ---> gamma-ray energy ---> final level energy ---> g.-ray intensity ---> g-ray int. uncertainty ---> level population int. ---> g.-ray Branching ratio
//see comments for ReadDecayScheme() for more details
vector<tuple<double, double, double, double, double, double >> AssignedTransition;


TGraph *gEff = new TGraph(); //graph for efficiency curve
TGraphErrors *gSigma = new TGraphErrors(); //graph for peak widths, !!!Please use sigma not FWHM
TGraphErrors *gEscPeaks = new TGraphErrors(); //graph to determine escape peak intensity

TF1 *fWidth; //linear function to fit peak widths as a function of energy
TF1 *fEscPeak;	//quadratic function to fit escape peak intensities as a function of energy

TFile *fRealData; //file contaaining real data
TH1D *hRealSpectra;	//real experimental histogram
TH1D *hSimPeaks; //simulated peak intensities
TH1D *hEscPeaks; //simulated escape peak intensities
TH1D *hBkgr; //simulated background extracted from real spectra
TH1D *hFullSim; //full simulated spectra including simulated full-energy and single-escape peaks on top of background 
TH1D *hSimPeaks_Bkgr;  //simulated spectra using only simulated full-energy peaks on top of background 
TH1D *hEscPeaks_Bkgr; //simulated spectra using only simulated single-escape peaks on top of background 

//this function calls for your input decay scheme
//This function requires a 5 column text file as the input
//The input is as follows
//Initial level energy ---> gamma-ray energy ---> final level energy ---> gamma-ray intensity ---> gamma-ray intensity uncertainty
//the gamma rays are added to the vector 'AssignedTransition'
//the 'level population int.' and 'g.-ray Branching ratio' components are only used for simulations of coincidence spectra
//these values are not read in by this function but are calculated by the function CalcLevelFeedingAndGammaBR()
//this function is included in the SimulateCoincidences.cc code

void ReadDecayScheme(string filename = "TransitionList.dat"){
	ifstream input( filename.c_str() );
	if( !input.is_open() ){
		cout << filename << " is not open!" << endl;
		return;
	}
	double a[5];
	input >> a[0] >> a[1] >> a[2] >> a[3] >> a[4];
	while( !input.eof() ){
		AssignedTransition.push_back(make_tuple(a[0], a[1], a[2], a[3], 0, 0 ));
		if(PrintReadData) cout << a[0] << "\t" <<  a[1] << "\t" <<  a[2] << "\t" <<  a[3] << "\t" <<  a[4] << endl;
		input >> a[0] >> a[1] >> a[2] >> a[3] >> a[4];
	}
}


//The GetEfficiency() function is used to read in your gamma-ray efficiency and add produce a graph of efficiency as a function of energy
//this function reads in a two column text file
//the first column is energy in keV
//the second column is the efficiency at that energy
//This efficiency graph is used to correct intensities when filling the simulated spectra
void GetEfficiency(string eff_filename = "MyExpEffnew.dat"){
	ifstream myfitresult( eff_filename.c_str() );
	if( !myfitresult.is_open() ){
		cout << eff_filename << " is not open!" << endl;
		return;
	}
	double a[2];
	myfitresult >> a[0] >> a[1];
	while( !myfitresult.eof() ){
		gEff->SetPoint( gEff->GetN(), a[0], a[1]);
		myfitresult >> a[0] >> a[1];
	}
}

//the GetPeakWidth() function is used to set the peak width as a function of energy for your simulated spectra
//This file requires a four column text file, the input should be the following
//gamma-ray energy ---> g.-ray energy error ---> peak width ---> peak width error
//please use term sigma! NOT FWHM!!!!
//this function fills a graph with peak widths as a function of energy
//a linear function is then fit to the data
//this function is used to get the peak widths for the simulated spectra
void GetPeakWidth(string peak_widths_filename = "PeakWidths.dat"){

	gSigma->SetName("gSigma");
	gSigma->SetMarkerStyle(20);
	gSigma->SetMarkerColor(kBlue);
	
	ifstream input( peak_widths_filename.c_str() );
	if( !input.is_open() ){
		cout << peak_widths_filename << " is not open!" << endl;
		return;
	}
	double a[4];
	input >> a[0] >> a[1] >> a[2] >> a[3];
	while( !input.eof() ){
		if( a[2] > a[3] ){
			gSigma->SetPoint( gSigma->GetN(), a[0], a[2]);
			gSigma->SetPointError( gSigma->GetN()-1,  a[1], a[3]);
			if(PrintPeakWidths) cout << a[0] << "\t" << a[1] << "\t" << a[2] << "\t" << a[3] << endl;
		}
			input >> a[0] >> a[1] >> a[2] >> a[3];
		
	}
	fWidth = new TF1("fWidth","[0]+[1]*x",0,8000);
	fWidth->SetParameters(9.57477e-01, 2.59267e-04);
	gSigma->Draw("AP");
	gSigma->Fit("fWidth");
	fWidth->Draw("same");
}

//the ReadEscapePeaks() function is used to get escape peak intensites relative to the full energy peak as a function of energy
//This file requires a four column text file, the input should be the following
//gamma-ray energy ---> g.-ray energy error ---> Esc.-Peak Int. / Full-Energy Peak Int  ---> uncertainty (Esc.-Peak Int. / Full-Energy Peak Int)
//this function is required for the simulation of escape peak intensities
void ReadEscapePeaks(string EscPeaksFilename = "EscapePeaks.dat"){

	gEscPeaks->SetName("gSigma");
	gEscPeaks->SetMarkerStyle(20);
	gEscPeaks->SetMarkerColor(kBlue);

	ifstream input( EscPeaksFilename.c_str() );
	if( !input.is_open() ){
		cout << EscPeaksFilename << " is not open!" << endl;
		return;
	}
	double a[4];
	input >> a[0] >> a[1] >> a[2] >> a[3];
	while( !input.eof() ){
		if( a[2] > a[3] ){
			gEscPeaks->SetPoint( gEscPeaks->GetN(), a[0], a[2]);
			gEscPeaks->SetPointError( gEscPeaks->GetN()-1,  a[1], a[3]);
			if(PrintEscPeakData) cout << a[0] << "\t" << a[1] << "\t" << a[2] << "\t" << a[3] << endl;
		}	
			input >> a[0] >> a[1] >> a[2] >> a[3];		
	}
	
	fEscPeak = new TF1("fEscPeak","[0]+[1]*x+[2]*x*x",0,8000);
	fEscPeak->SetParameters(-1.45859e-02, 1.11648e-06, 7.51546e-09);
	gEscPeaks->Draw("AP");
	gEscPeaks->Fit("fEscPeak");
	fEscPeak->Draw("same");
}

//this function is used to extract a background from your experimntal spectrum
//the default values of the function will extract a rather crude background
//the user is free to play with these parameters to try and improve the extracted background
void GetSpectrumBackground(TH1D *h, int iterations = 50, int decreasewindow = 1, int backorder = 2, bool smoothing = false, int smoothwindow = 3, bool compton = 0){

	//extract binning information from the experimental spectrum
	//we want to use the same binning to compare
	const int Nbins = h->GetXaxis()->GetNbins();
	double x_low = h->GetXaxis()->GetBinLowEdge(1);
	double x_max = h->GetXaxis()->GetBinUpEdge(Nbins);
	double source[Nbins];
	
	hBkgr = new TH1D("hBkgr","Simulated Background",Nbins,x_low,x_max);
	
	TSpectrum *s = new TSpectrum();
	for (int i = 0; i < Nbins; i++) source[i]=h->GetBinContent(i + 1);
	//s->Background(source,Nbins,75,TSpectrum::kBackDecreasingWindow, TSpectrum::kBackOrder2,kFALSE, TSpectrum::kBackSmoothing3,kFALSE);
	s->Background(source,Nbins,iterations,decreasewindow,backorder,smoothing,smoothwindow,compton);
	for (int i = 0; i < Nbins; i++) hBkgr->SetBinContent(i + 1,source[i]);      
	h->Draw("hist");
	hBkgr->SetLineColor(kOrange+1);
	hBkgr->Draw("SAME L");
                 
}

//this function is used to get your real experimental spectra
//user should provide the name of the root file and the name of the spectrum
void GetRealSpectra(string rootfilename = "ExampleFile.root", string RealHistName = "hgE_56Co"){

	//open root file
	if(gSystem->AccessPathName(rootfilename.c_str())){ //check to see root file exists!
		cout << "The file " << rootfilename << " doesn't exists\n";
		return;
	}
	else fRealData = TFile::Open( rootfilename.c_str() ); //opening root file
	
	if( !fRealData->GetListOfKeys()->Contains( RealHistName.c_str() )){ //check to see histogram exists!
		cout << "The histogram " << RealHistName << " doesn't exists\n";
		return;
	}	
	else{ //getting root spectrum
		hRealSpectra = (TH1D*)fRealData->Get( RealHistName.c_str() );
		hRealSpectra->Draw("hist");
		hRealSpectra->GetXaxis()->SetTitle("Energy (keV)");
		hRealSpectra->GetXaxis()->CenterTitle();
	}
}

//vector<tuple<double, double, double, double>> EscPeaksData;

//Addition info for GetSpectrumBackground() function
/*
       kBackOrder2 =0,
       kBackOrder4 =1,
       kBackOrder6 =2,
       kBackOrder8 =3,
       kBackIncreasingWindow =0,
       kBackDecreasingWindow =1,
       kBackSmoothing3 =3,
       kBackSmoothing5 =5,
       kBackSmoothing7 =7,
       kBackSmoothing9 =9,
       kBackSmoothing11 =11,
       kBackSmoothing13 =13,
       kBackSmoothing15 =15
       
       Background(Double_t* spectrum, Int_t ssize, Int_t numberIterations, Int_t direction, Int_t filterOrder, bool smoothing, Int_t smoothWindow, bool compton)
*/
