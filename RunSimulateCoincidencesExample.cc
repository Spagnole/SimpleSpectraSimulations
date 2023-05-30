#include "SimulateCoincidences.cc"


void RunSim_Gate1238(){

	ReadDecayScheme("TransitionList.dat");
	GetEfficiency("MyExpEffnew.dat");
	GetPeakWidth("PeakWidths.dat");
	ReadEscapePeaks("EscapePeaks.dat");
	GetRealSpectra("ExampleFile.root","hggE_Gate1238");	
	GetSpectrumBackground(hRealSpectra,50);
	
	GetLevelList();
	FixTransitionList();
	CalcLevelFeedingAndGammaBR();
	FindCoincidences();
	ReduceGammaGammaList();
	FillCoincMatrix();
	
	GateOnSimMat(1238-5,1238+5);
	FillEscPeakSpec(SimMatProj);
	BuildSimuledSpectra(SimMatProj,hBkgr,hEscPeak);
	
}

void RunSim_Gate1038(){

	ReadDecayScheme("TransitionList.dat");
	GetEfficiency("MyExpEffnew.dat");
	GetPeakWidth("PeakWidths.dat");
	ReadEscapePeaks("EscapePeaks.dat");
	GetRealSpectra("ExampleFile.root","hggE_Gate1038");	
	GetSpectrumBackground(hRealSpectra,50);
	
	GetLevelList();
	FixTransitionList();
	CalcLevelFeedingAndGammaBR();
	FindCoincidences();
	ReduceGammaGammaList();
	FillCoincMatrix();
	
	GateOnSimMat(1038-5,1038+5);
	FillEscPeakSpec(SimMatProj);
	BuildSimuledSpectra(SimMatProj,hBkgr,hEscPeak);
	
}

void RunSim_Gate2034(){

	ReadDecayScheme("TransitionList.dat");
	GetEfficiency("MyExpEffnew.dat");
	GetPeakWidth("PeakWidths.dat");
	ReadEscapePeaks("EscapePeaks.dat");
	GetRealSpectra("ExampleFile.root","hggE_Gate2034");	
	GetSpectrumBackground(hRealSpectra,50);
	
	GetLevelList();
	FixTransitionList();
	CalcLevelFeedingAndGammaBR();
	FindCoincidences();
	ReduceGammaGammaList();
	FillCoincMatrix();
	
	GateOnSimMat(2034-5,2034+5);
	FillEscPeakSpec(SimMatProj);
	BuildSimuledSpectra(SimMatProj,hBkgr,hEscPeak);
	
}

