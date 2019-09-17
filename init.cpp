
void init(){

	gROOT->ProcessLine(".L Risultati.h+");
	gROOT->ProcessLine(".L Risultati.cpp+");
	gROOT->ProcessLine(".L Header.h+");
	gROOT->ProcessLine(".L Header.cpp+");
	gROOT->ProcessLine(".L istogramma.cpp");
	gROOT->ProcessLine(".L fit.cpp");
	gROOT->ProcessLine(".L fft.cpp");
	gROOT->ProcessLine(".L LetturaFile.cpp");
	gROOT->ProcessLine(".L GetHeader.cpp");
	gROOT->ProcessLine(".L txtToROOT.cpp");
	gROOT->ProcessLine(".L readcsv.C");
	gROOT->ProcessLine(".L calibrazione.cpp");

}
