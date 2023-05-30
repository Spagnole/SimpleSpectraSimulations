#include "SimulateSingles.cc"

void RunFullSimulation(){

	ReadDecayScheme("TransitionList.dat");
	GetEfficiency("MyExpEffnew.dat");
	GetPeakWidth("PeakWidths.dat");
	ReadEscapePeaks("EscapePeaks.dat");
	GetRealSpectra("ExampleFile.root","hgE_56Co");
	GetSpectrumBackground(hRealSpectra);
	FillSimulation();
	FillEscapePeaks();
	BuildSimuledSpectra(hSimPeaks,hBkgr,hEscPeaks);
	
}
