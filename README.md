# adc-test

Steps:
	1) run ".x init.cpp" in the ROOT shell in order to load the relevant files;
	2) use the testcsv() function in readcsv.C to convert csv to txt; 
	3) use txtToROOT() in txtToROOT.cpp to convert the txt file to a ROOT file with a TGraph containing the samples and a TH1D with the codes distribution;
	4) (optional) run calibrazione() to quantify and correct interleaving errors. The calibration parameters can be saved in the ROOT file;
	5) run the relevant analysis with fit(), fft() or istogramma().

Note on coding style: the indentation character is tab.
