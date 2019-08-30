
#include "Risultati.h"
#include "Header.h"

#if !defined (__CLING__) || defined (__ROOTCLING__)
//		C++
#include <cmath>
#include <iostream>
#include <unordered_map>
//		ROOT
#include "TH1F.h"
#include "TF1.h"
#include "TMath.h"
#include <TCanvas.h>
#include <TGraph.h>
#include <TPad.h>
#include <TStyle.h>
#include <TH2D.h>
#endif

//NOTA BENE: se non si conosce bene l'ampiezza dell'onda allora l'SNDR è inaccurato!

//	<DA FARSI>
//	</DA FARSI>

int fit(

	std::string nomeFile = "datiFFiT.root",
	bool salva = 0,
	double MSPS = 160,
	bool useCalibratedSamples = 1,		//analizza i dati calibrati, se presenti
	int split = 0 		//0: tutti, 1: dispari, 2: pari

	){

	//definizioni iniziali parametri ADC
	const double nBit = 12;
	const int nCh = (int) TMath::Power(2, nBit);
	const int chMax = nCh - 1;
	const double vmax = 0.6;
	const double vmin = -0.6;
	const double FSR = vmax - vmin;
	const double LSB = FSR / nCh;
	//~ const double MSPS = 160;
	if ( split != 0 ) MSPS = MSPS / 2;
	const double fSamp = MSPS * 1000000;
	const double tSamp = 1 / fSamp;

	//definizioni varie
	//~ std::string nomeTGraph = "tGr";
	const std::string nomeTGraphDefault = "grSamp";
	std::string nomeTGraph = nomeTGraphDefault;
	if ( split == 1 && useCalibratedSamples == 0 ) nomeTGraph += "Odd";
	if ( split == 2 && useCalibratedSamples == 0 ) nomeTGraph += "Even";
	if ( useCalibratedSamples ) nomeTGraph += "Calib";
	const std::string ampKey = "amp[V]";
	const std::string freqKey = "freq[Hz]";
	std::string sampleType = "";
	if ( !useCalibratedSamples && !split ) sampleType += "All";
	if ( split == 1 ) sampleType += "Odd";
	if ( split == 2 ) sampleType += "Even";
	if ( useCalibratedSamples ) sampleType += "Calibrated";

	const bool grafica = 1;
	const bool output = 1;

	//	<LETTURA FILE>

	TFile tf1( nomeFile.c_str(), "update");

	//controllo esistenza del file
	if ( !tf1.IsOpen() ){

		std::cout << "File \"" << nomeFile.c_str() << "\" not found" << std::endl;
		std::cout << "Aborting execution" << std::endl;
		return 1;

	}

	TGraph * gr1 = (TGraph *) tf1.Get( nomeTGraph.c_str() );

	if ( gr1 == nullptr && useCalibratedSamples ) {

		std::cout << "TGraph \"" << nomeTGraph.c_str() << "\" not found " << std::endl;
		std::cout << "Reverting to \"" << nomeTGraphDefault.c_str() << "\" " << std::endl;
		nomeTGraph = nomeTGraphDefault;
		gr1 = (TGraph *) tf1.Get( nomeTGraph.c_str() );

	}

	if ( gr1 == nullptr ){

		std::cout << "TGraph \"" << nomeTGraph.c_str() << "\" not found " << std::endl;
		std::cout << "Aborting execution" << std::endl;
		return 2;

	}

	gr1->SetMarkerStyle(7);

	//		<LETTURA HEADER>
	Header * headerObj = (Header *) tf1.Get( "Header" );

	if ( headerObj == nullptr ){

		std::cout << "Header not found" << std::endl;
		std::cout << "Aborting execution" << std::endl;
		delete gr1;
		tf1.Close();
		return 3;

	}

	std::unordered_map < std::string, double > headerMap = headerObj->GetHeader();
	delete headerObj;

	double amp = 0;
	double freq = 0;

	try {

		amp = headerMap.at( ampKey );
		freq = headerMap.at( freqKey );

	} catch ( const std::out_of_range& invArg ) {

		std::cout << "Key not found" << std::endl;
		std::cout << "Aborting execution" << std::endl;
		delete gr1;
		tf1.Close();
		return 4;

	}

	//		</LETTURA HEADER>

	//altre definizioni
	double omega = 2 * TMath::Pi() * freq;
	double T = 1 / freq;
	double Tsamp = T / tSamp;	//periodo in unità di campionamenti ideali
	double omegaSamp = 2 * TMath::Pi() / Tsamp;

	int nPts = gr1->GetN();

	double * iSamp = new double[ nPts ];
	double * res = new double[ nPts];
	double * ch = new double[ nPts ];

	//	</LETTURA FILE>

	//rette
	TF1	* line1 = new TF1( "line1", "0.5", 0, nPts - 1 );
	line1->SetLineColor(kOrange - 3);
	line1->SetLineWidth(1);
	line1->SetLineStyle(1);

	TF1	* line2 = new TF1( "line2", "-0.5", 0, nPts - 1 );
	line2->SetLineColor(kOrange - 3);
	line2->SetLineWidth(1);
	line2->SetLineStyle(1);

	TF1	* line3 = new TF1( "line3", "0", 0, nPts - 1 );
	line3->SetLineColor(kOrange - 3);
	line3->SetLineWidth(1);
	line3->SetLineStyle(1);

	TF1	* line4 = new TF1( "line4", std::to_string(nCh).c_str(), 0, nPts - 1 );
	line4->SetLineColor(kOrange - 3);
	line4->SetLineWidth(1);
	line4->SetLineStyle(1);

	//istogramma distribuzione residui
	TH1F histRes( "histRes", "histRes", 800, -20, 20);

	//istogramma campionamenti in funzione della fase
	TH1D histPhase( "histPhase", "Sampling phase", nPts, 0 - 2 * TMath::Pi() / ( 2 * nPts ) , 2 * TMath::Pi() + 2 * TMath::Pi() / ( 2 * nPts ));

	//scatter plot dei residui in funzione della fase
	TH2D histErrPhase("histErrPhase", "Phase-residuals distribution", 50, 0, 2*TMath::Pi(), 5000, -50, 50 );

	//ricava valori iniziali dei parametri del fit dai dati
	int nMin_exp = TMath::LocMin( nPts, gr1->GetY() );
	int nMax_exp = TMath::LocMax( nPts, gr1->GetY() );
	double chMin_exp = gr1->GetY()[ nMin_exp ];
	double chMax_exp = gr1->GetY()[ nMax_exp ];
	double A_exp = ( chMax_exp - chMin_exp ) / 2;
	double offset_exp = chMax_exp - A_exp;

	//fit
	TF1 * sinFit = new TF1( "sinFit", "[0] * TMath::Sin( [1] * x + [2] ) + [3]", 0, nPts - 1 );
	sinFit->SetParameter(0, A_exp);
	sinFit->SetParLimits(0, A_exp*0.95, A_exp*1.05);
	sinFit->SetParName(0, "Amplitude");
	sinFit->FixParameter( 1, omegaSamp );		//importante, sennò non riesce a fare il fit
	//~ sinFit->SetParameter( 1, omegaSamp );		//importante, sennò non riesce a fare il fit
	//~ sinFit->SetParLimits(1, omegaSamp * 0.999, omegaSamp * 1.001);
	sinFit->SetParName(1, "Omega");
	//~ sinFit->SetParameter(2, phi);
	sinFit->SetParName(2, "Phase");
	sinFit->SetParameter(3, offset_exp);
	sinFit->SetParLimits(3, offset_exp * 0.95, offset_exp * 1.05);
	sinFit->SetParName(3, "Offset");
	//~ gr1->Fit("sinFit", "RN0Q");		//fithttps://stackoverflow.com/questions/160930/how-do-i-check-if-an-integer-is-even-or-odd
	gr1->Fit("sinFit", "N0QB");		//fit

	//calcolo varianza
	double var = 0;

	for ( long i = 0; i < nPts; i++ ){

		double x, y = 0;
		gr1->GetPoint(i, x, y);
		iSamp[i] = x;
		res[i] = y - sinFit->Eval( x );
		var += res[i] * res[i];
		histRes.Fill( res[i] );

		//fase
		double fase = fmod( sinFit->GetParameter(1) * i + sinFit->GetParameter(2),  2 * TMath::Pi() );
		histPhase.Fill( fase );

		//scatter plot
		histErrPhase.Fill(fase, res[i]);

	}

	//grafico dei residui in funzione del tempo
	TGraph gr2( nPts, iSamp, res );
	gr2.SetName( "grResT" );

	//grafico dei residui tra fit e valore vero
	// TGraph gr3( nPts, t, res2);
	// gr3.SetMarkerStyle(6);
	// gr3.SetLineStyle(3);

	//varianze
	var = var * LSB * LSB / nPts;
	double varQ = ( LSB * LSB ) / 12;
	double sigmaQ = TMath::Sqrt( varQ );
	double NAD = TMath::Sqrt(var);

	//ENOB:
	double enob1 = nBit - TMath::Log2( NAD / sigmaQ );
	double enob2 = TMath::Log2(  FSR / ( TMath::Sqrt(12) * NAD )  );
	//~ double enob3 = nBit - TMath::Log2( histRes.GetRMS()*LSB / sigmaQ );

	//SNDR:
	double Arms = amp / TMath::Sqrt(2);
	double sinad1 = 20 * TMath::Log10( Arms / NAD );
	//~ double sinad2 = (enob1 * 6.02 + 1.76);

	if (output){

		std::cout << "======================================= Results =======================================" << std::endl;

		std::cout << "Dataset: " << nomeFile.c_str() << std::endl;
		std::cout << "Type: " << sampleType.c_str() << std::endl;

		std::cout << "Sine fit parameters:" << std::endl;
		std::cout << "\t0) Amplitude: " << sinFit->GetParameter(0) << std::endl;
		std::cout << "\t1) Frequency [Hz]: " << sinFit->GetParameter(1) / ( 2 * TMath::Pi() ) * fSamp << std::endl;	//è in unità di campionamenti
		std::cout << "\t2) Phase: " << sinFit->GetParameter(2) << std::endl;
		std::cout << "\t3) Offset: " << sinFit->GetParameter(3) << std::endl;

		std::cout << "Variance [LSB]: " << var / ( LSB * LSB ) << std::endl;
		std::cout << "Variance (ideal) [LSB]: " << varQ / ( LSB * LSB ) << std::endl;

		//~ std::cout << "ENOB, metodo 1: " << enob1 << std::endl;
		std::cout << "ENOB [bit]: " << enob2 << std::endl;
		//~ // std::cout << "ENOB, metodo 3: " << enob3 << std::endl;

		std::cout << "SNDR [dB]: " << sinad1 << std::endl;	//questo è il metodo presente nello standard IEEE 2017
		//~ std::cout << "SINAD[dB], direttamente da ENOB:" << sinad2 << std::endl;

		std::cout << "========================================================================================" << std::endl;

	}

	// <OUTPUT>

	if (grafica || salva){

		if (split == 0) gr1->SetTitle("Sampled signal");
		if (split == 1) gr1->SetTitle("Sampled signal, odd samples");
		if (split == 2) gr1->SetTitle("Sampled signal, even samples");
		if (split == 0) gr1->GetXaxis()->SetTitle("Sample");
		if (split == 1) gr1->GetXaxis()->SetTitle("Sample (odd)");
		if (split == 2) gr1->GetXaxis()->SetTitle("Sample (even)");
		gr1->GetYaxis()->SetTitle("ADC code");
		gr1->GetXaxis()->SetRangeUser( 0, nPts * 0.05);
		gr1->GetYaxis()->SetRangeUser( -100, nCh + 100);

		gr2.SetTitle("Residuals as a function of time");
		if (split == 0) gr2.GetXaxis()->SetTitle("Sample");
		if (split == 1) gr2.GetXaxis()->SetTitle("Sample (odd)");
		if (split == 2) gr2.GetXaxis()->SetTitle("Sample (even)");
		gr2.GetYaxis()->SetTitle("#Delta[ADC counts]");
		gr2.GetXaxis()->SetRangeUser( -1, nPts );
		gr2.SetLineStyle(1);
		gr2.SetLineWidth(1);
		gr2.SetMarkerStyle(6);
		// gr2.SetMarkerStyle(7);
		// gr2.SetMarkerSize(0.5);

		histRes.SetTitle( "Distribution of residuals" );
		histRes.GetXaxis()->SetTitle( "#Delta[ADC counts]" );
		histRes.GetYaxis()->SetTitle( "Counts" );
		histRes.GetXaxis()->SetRangeUser( histRes.GetMean() - 5*histRes.GetStdDev(), histRes.GetMean() + 5*histRes.GetStdDev() );
		//~ histRes.SetLineColor( kRed + 2 );
		//~ histRes.SetFillColor( kOrange - 3 );
		//~ histRes.SetFillStyle( 3004 );
		//~ histRes.GetXaxis()->SetNdivisions(20);

		histPhase.SetLineColor( kRed + 2 );
		//~ histPhase.SetFillColor( kOrange - 3 );
		//~ histPhase.SetFillStyle( 3004 );
		histPhase.GetXaxis()->SetTitle("Phase[rad]");
		histPhase.GetYaxis()->SetTitle("Counts");
		histPhase.SetMarkerStyle(7);

		histErrPhase.GetXaxis()->SetTitle("Phase[rad]");
		histErrPhase.GetYaxis()->SetTitle("#Delta[ADC counts]");
		histErrPhase.GetYaxis()->SetRangeUser( histRes.GetMean() - 5*histRes.GetStdDev(), histRes.GetMean() + 5*histRes.GetStdDev() );

		sinFit->SetLineWidth(1);
		sinFit->SetLineColor( kRed );

	}

	if (grafica){

		TCanvas * c1 = new TCanvas("c1", "", 1250, 650);
		c1->Divide(2, 3);

		//grafico campionamenti
		c1->cd(1);
		gPad->SetGrid();
		gr1->DrawClone("AP");
		line3->DrawCopy("same");
		line4->DrawCopy("same");

		//grafico residui
		c1->cd(3);
		gPad->SetGrid();
		gr2.DrawClone("LPA");
		line1->DrawCopy("same");
		line2->DrawCopy("same");

		//distribuzione dei residui
		c1->cd(4);
		gPad->SetGrid();
		histRes.DrawCopy();

		c1->Update();

		c1->cd(5);
		gPad->SetGrid();
		gStyle->SetOptStat(0);
		histPhase.DrawCopy("");

		c1->Update();

		//scatter plot
		c1->cd(6);
		gStyle->SetOptStat(0);
		gPad->SetGrid();
		histErrPhase.DrawCopy("COLZ");

		c1->Update();

		//grafico fit; deve essere l'ultimo!
		gr1->SetTitle("Sinewave fit");
		c1->cd(2);
		gPad->SetGrid();
		gStyle->SetOptStat(1110);
		gStyle->SetOptFit(1110);
		gr1->DrawClone("PA");
		sinFit->DrawCopy("same");
		line3->DrawCopy("same");
		line4->DrawCopy("same");
		gPad->Modified();

	}

	// <SALVATAGGIO>

	if (salva){

		Risultati * results = (Risultati *) tf1.Get("Risultati");
		if ( results == nullptr ) results = new Risultati();
		results->AddFitData( "Variance", var / ( LSB * LSB ) );
		results->AddFitData( "ENOB", enob2 );
		results->AddFitData( "SNDR", sinad1 );
		results->Write("", TObject::kOverwrite);
		delete results;

		gr2.Write("", TObject::kOverwrite);
		sinFit->Write("", TObject::kOverwrite);
		histRes.Write("", TObject::kOverwrite);
		histPhase.Write("", TObject::kOverwrite);
		histErrPhase.Write("", TObject::kOverwrite);

	}

	// </SALVATAGGIO>

	// </OUTPUT>

	delete sinFit;
	delete line1, line2, line3, line4;
	delete gr1;
	delete [] iSamp;
	delete [] ch;
	delete [] res;

	tf1.Close();

	return 0;

}
