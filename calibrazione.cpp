
#include "Risultati.h"
#include "Header.h"

#if !defined (__CINT__) || defined (__MAKECINT__)
#include "TH1F.h"
#include "TF1.h"
#include "TMath.h"
#include <math.h>
#include <iostream>
#include <cmath>	//per ceil
#include <unordered_map>
#include <TCanvas.h>
#include <TGraph.h>
#include <TPad.h>
#include <TStyle.h>
#include <TH2D.h>
#endif

//	<DA FARSI>
//	-aggiungere calibrazione di fase
//	</DA FARSI>

int calibrazione(

	std::string nomeFile = "datiFFiT.root",
	bool salva = 0,
	bool grafica = 1,
	bool output = 1,
	double MSPS = 160

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
	const double fSamp = MSPS * 1000000;
	const double tSamp = 1 / fSamp;

	//definizioni varie
	//~ std::string nomeTGraph = "tGr";
	const std::string nomeTGraph = "grSamp";
	const std::string nomeTGraphCalib = nomeTGraph + "Calib";
	const std::string ampKey = "amp[V]";
	const std::string freqKey = "freq[Hz]";

	//	<LETTURA FILE>

	TFile tf1( nomeFile.c_str(), "update");

	//controllo esistenza del file
	if ( !tf1.IsOpen() ){

		return 1;

	}

	TGraph * gr1 = (TGraph *) tf1.Get( nomeTGraph.c_str() );

	if ( gr1 == nullptr ){

		std::cout << "TGraph \"" << nomeTGraph.c_str() << "\" non trovato " << std::endl;
		return 2;

	}

	gr1->SetMarkerStyle(7);

	//		<LETTURA HEADER>
	Header * headerObj = (Header *) tf1.Get( "Header" );

	if ( headerObj == nullptr ){

		std::cout << "Header non trovato" << std::endl;
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

		std::cout << "Chiave non trovata" << std::endl;
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
	TF1 * line1 = new TF1( "line1", "0.5", 0, nPts - 1 );
	line1->SetLineColor(kOrange - 3);
	line1->SetLineWidth(1);
	line1->SetLineStyle(1);

	TF1 * line2 = new TF1( "line2", "-0.5", 0, nPts - 1 );
	line2->SetLineColor(kOrange - 3);
	line2->SetLineWidth(1);
	line2->SetLineStyle(1);

	TF1 * line3 = new TF1( "line3", "0", 0, nPts - 1 );
	line3->SetLineColor(kOrange - 3);
	line3->SetLineWidth(1);
	line3->SetLineStyle(1);

	TF1 * line4 = new TF1( "line4", std::to_string(nCh).c_str(), 0, nPts - 1 );
	line4->SetLineColor(kOrange - 3);
	line4->SetLineWidth(1);
	line4->SetLineStyle(1);

	//<FIT INTERLEAVING>

	TGraph * grEven = new TGraph();
	TGraph * grOdd = new TGraph();

	grEven->SetTitle("Even samples");
	grEven->GetXaxis()->SetTitle("Sample number");
	grEven->GetYaxis()->SetTitle("ADC count");
	grOdd->SetTitle("Odd samples");
	grOdd->GetXaxis()->SetTitle("Sample number");
	grOdd->GetYaxis()->SetTitle("ADC count");

	int nEven = 0;
	int nOdd = 0;

	//riempimento
	for ( long i = 0; i < nPts; i++ ) {

		int yTemp = gr1->GetY()[i];

		if ( i % 2 == 0 ) {

			grEven->SetPoint(nEven, i, yTemp);
			nEven++;

		} else {

			grOdd->SetPoint(nOdd, i, yTemp);
			nOdd++;

		}

	}

	//even
	int nMinEven_exp = TMath::LocMin( nEven, grEven->GetY() );
	int nMaxEven_exp = TMath::LocMax( nEven, grEven->GetY() );
	double chMinEven_exp = grEven->GetY()[ nMinEven_exp ];
	double chMaxEven_exp = grEven->GetY()[ nMaxEven_exp ];
	double AEven_exp = ( chMaxEven_exp - chMinEven_exp ) / 2;
	double offsetEven_exp = chMaxEven_exp - AEven_exp;

	//odd
	int nMinOdd_exp = TMath::LocMin( nOdd, grOdd->GetY() );
	int nMaxOdd_exp = TMath::LocMax( nOdd, grOdd->GetY() );
	double chMinOdd_exp = grOdd->GetY()[ nMinOdd_exp ];
	double chMaxOdd_exp = grOdd->GetY()[ nMaxOdd_exp ];
	double AOdd_exp = ( chMaxOdd_exp - chMinOdd_exp ) / 2;
	double offsetOdd_exp = chMaxOdd_exp - AOdd_exp;

	//even
	TF1 * sinEven = new TF1( "sinEven", "[0] * TMath::Sin( [1] * x + [2] ) + [3]", 0, nPts - 1 );
	sinEven->SetParameter(0, AEven_exp);
	sinEven->SetParLimits(0, AEven_exp*0.95, AEven_exp*1.05);
	sinEven->SetParName(0, "Amplitude");
	sinEven->FixParameter( 1, omegaSamp );		//importante, sennò non riesce a fare il fit
	//~ sinEven->SetParLimits(1, omegaSamp * 0.95, omegaSamp * 1.05);
	sinEven->SetParName(1, "Omega");
	sinEven->SetParName(2, "Phase");
	sinEven->SetParameter(3, offsetEven_exp);
	sinEven->SetParLimits(3, offsetEven_exp * 0.95, offsetEven_exp * 1.05);
	sinEven->SetParName(3, "Offset");
	grEven->Fit("sinEven", "N0QB");		//fit

	//odd
	TF1 * sinOdd = new TF1( "sinOdd", "[0] * TMath::Sin( [1] * x + [2] ) + [3]", 0, nPts - 1 );
	sinOdd->SetParameter(0, AOdd_exp);
	sinOdd->SetParLimits(0, AOdd_exp*0.95, AOdd_exp*1.05);
	sinOdd->SetParName(0, "Amplitude");
	sinOdd->FixParameter( 1, omegaSamp );		//importante, sennò non riesce a fare il fit
	//~ sinOdd->SetParLimits(1, omegaSamp * 0.95, omegaSamp * 1.05);
	sinOdd->SetParName(1, "Omega");
	sinOdd->SetParName(2, "Phase");
	sinOdd->SetParameter(3, offsetOdd_exp);
	sinOdd->SetParLimits(3, offsetOdd_exp * 0.95, offsetOdd_exp * 1.05);
	sinOdd->SetParName(3, "Offset");
	grOdd->Fit("sinOdd", "N0QB");		//fit

	//</FIT INTERLEAVING>

	//<CALIBRAZIONE>
	TGraph * grCalib = new TGraph();
	grCalib->SetName( nomeTGraphCalib.c_str() );

	//fattori di calibrazione
	double deltaOffset = 0;
	double amplitudeRatio = 0;
	bool calibrationParity = 0;		//0 pari, 1 dispari
	int halfChannel = round ( nCh / 2 );

	double AOdd = sinOdd->GetParameter(0);
	double AEven = sinEven->GetParameter(0);
	double offsetOdd = sinOdd->GetParameter(3);
	double offsetEven = sinEven->GetParameter(3);

	//scelta su quale core effettuare la calibrazione
	if (  TMath::Abs( offsetEven - halfChannel ) > TMath::Abs( offsetOdd - halfChannel )  ) calibrationParity = 0;
	if ( AOdd > AEven ) calibrationParity = 0;

	for ( int i = 0; i < nPts; i++ ) {

		bool sampleParity = ( i % 2 == 0 ) ? 0 : 1;		//0 pari, 1 dispari

		if ( sampleParity == calibrationParity ){

			double uncalibADC = gr1->GetY()[i];
			int calibADC = 0;

			if ( calibrationParity == 0 ) {

				uncalibADC = uncalibADC - offsetEven;
				uncalibADC = uncalibADC * ( AOdd / AEven );
				uncalibADC = uncalibADC + offsetOdd;

			} else {

				uncalibADC = uncalibADC - offsetOdd;
				uncalibADC = uncalibADC * ( AEven / AOdd );
				uncalibADC = uncalibADC + offsetEven;

			}

			calibADC = uncalibADC;
			grCalib->SetPoint(i, i, calibADC);

		} else {

			grCalib->SetPoint(i, i, gr1->GetY()[i]);

		}

	}

	//</CALIBRAZIONE>

	//ricava valori iniziali dei parametri del fit dai dati
	int nMin_exp = TMath::LocMin( nPts, gr1->GetY() );
	int nMax_exp = TMath::LocMax( nPts, gr1->GetY() );
	double chMin_exp = gr1->GetY()[ nMin_exp ];
	double chMax_exp = gr1->GetY()[ nMax_exp ];
	double A_exp = ( chMax_exp - chMin_exp ) / 2;
	double offset_exp = chMax_exp - A_exp;

	//fit
	TF1 * sinFitCalib = new TF1( "sinFitCalib", "[0] * TMath::Sin( [1] * x + [2] ) + [3]", 0, nPts - 1 );
	sinFitCalib->SetParameter(0, A_exp);
	sinFitCalib->SetParLimits(0, A_exp*0.95, A_exp*1.05);
	sinFitCalib->SetParName(0, "Amplitude");
	sinFitCalib->FixParameter( 1, omegaSamp );		//importante, sennò non riesce a fare il fit
	//~ sinFitCalib->SetParLimits(1, omegaSamp * 0.95, omegaSamp * 1.05);
	sinFitCalib->SetParName(1, "Omega");
	//~ sinFitCalib->SetParameter(2, phi);
	sinFitCalib->SetParName(2, "Phase");
	sinFitCalib->SetParameter(3, offset_exp);
	sinFitCalib->SetParLimits(3, offset_exp * 0.95, offset_exp * 1.05);
	sinFitCalib->SetParName(3, "Offset");
	//~ gr1->Fit("sinFitCalib", "RN0Q");		//fit
	grCalib->Fit("sinFitCalib", "N0QB");		//fit

	if (output){

		std::cout << "===================== Results =====================" << std::endl;

		std::cout << "Even samples, sine fit parameters:" << std::endl;
		std::cout << "\t0) Amplitude: " << sinEven->GetParameter(0) << std::endl;
		std::cout << "\t1) Frequency [Hz]: " << sinEven->GetParameter(1) / ( 2 * TMath::Pi() ) * fSamp << std::endl;	//è in unità di campionamenti
		std::cout << "\t2) Phase: " << sinEven->GetParameter(2) << std::endl;
		std::cout << "\t3) Offset: " << sinEven->GetParameter(3) << std::endl;

		std::cout << "Odd samples, sine fit parameters:" << std::endl;
		std::cout << "\t0) Amplitude: " << sinOdd->GetParameter(0) << std::endl;
		std::cout << "\t1) Frequency [Hz]: " << sinOdd->GetParameter(1) / ( 2 * TMath::Pi() ) * fSamp << std::endl;	//è in unità di campionamenti
		std::cout << "\t2) Phase: " << sinOdd->GetParameter(2) << std::endl;
		std::cout << "\t3) Offset: " << sinOdd->GetParameter(3) << std::endl;

		std::cout << "Sine fit parameters after calibration:" << std::endl;
		std::cout << "\t0) Amplitude: " << sinFitCalib->GetParameter(0) << std::endl;
		std::cout << "\t1) Frequency [Hz]: " << sinFitCalib->GetParameter(1) / ( 2 * TMath::Pi() ) * fSamp << std::endl;	//è in unità di campionamenti
		std::cout << "\t2) Phase: " << sinFitCalib->GetParameter(2) << std::endl;
		std::cout << "\t3) Offset: " << sinFitCalib->GetParameter(3) << std::endl;

		std::cout << "Calibration parameters:" << std::endl;
		std::cout << "\tOffset difference (odd - even): " << offsetOdd - offsetEven << " ch" << std::endl;
		std::cout << "\tAmplitude difference (odd - even): " << AOdd - AEven << " ch" << std::endl;
		std::cout << "\tPhase difference (odd - even): " << sinOdd->GetParameter(2) - sinEven->GetParameter(2) << " rad" << std::endl;

		std::cout << "====================================================" << std::endl;

	}

	// <OUTPUT>

	if (grafica || salva){

		gr1->SetTitle("Uncalibrated samples");
		gr1->GetXaxis()->SetTitle("Sample");
		gr1->GetYaxis()->SetTitle("ADC count");
		gr1->GetXaxis()->SetRangeUser( 0, 0.005 * nPts );
		gr1->GetYaxis()->SetRangeUser( -100, nCh + 100);
		gr1->SetMarkerStyle(25);
		gr1->SetMarkerSize(0.5);
		//~ gr1->SetMarkerColor(kRed);

		grCalib->SetTitle("Calibrated samples");
		grCalib->GetXaxis()->SetTitle("Sample");
		grCalib->GetYaxis()->SetTitle("ADC count");
		grCalib->GetXaxis()->SetRangeUser( 0, 0.005 * nPts );
		grCalib->GetYaxis()->SetRangeUser( -100, nCh + 100);
		grCalib->SetMarkerStyle(24);
		grCalib->SetMarkerSize(0.5);
		//~ grCalib->SetMarkerColor(kGreen);

		grEven->SetMarkerStyle(25);
		grEven->SetMarkerSize(0.5);

		grOdd->SetMarkerStyle(25);
		grOdd->SetMarkerSize(0.5);

		sinFitCalib->SetLineWidth(1);
		sinFitCalib->SetLineColor( kRed );

		sinEven->SetLineWidth(1);
		//~ sinEven->SetLineStyle(3);
		sinEven->SetLineColor(kOrange);

		sinOdd->SetLineWidth(1);
		//~ sinOdd->SetLineStyle(3);
		sinOdd->SetLineColor(kGreen);

	}

	if (grafica){

		TCanvas * c1 = new TCanvas("c1", "", 1250, 650);
		c1->Divide(1, 2);

		c1->cd(1);
		gPad->SetGrid();
		gr1->DrawClone("AP");
		sinEven->DrawCopy("same");
		sinOdd->DrawCopy("same");
		line3->DrawCopy("same");
		line4->DrawCopy("same");

		c1->cd(2);
		gPad->SetGrid();
		//~ gr1->DrawClone("AP");
		grCalib->DrawClone("AP");
		sinFitCalib->DrawCopy("same");
		sinEven->DrawCopy("same");
		sinOdd->DrawCopy("same");
		line3->DrawCopy("same");
		line4->DrawCopy("same");

	}

	// <SALVATAGGIO>

	if (salva){

		grCalib->Write("", TObject::kOverwrite);
		sinEven->Write("", TObject::kOverwrite);
		sinOdd->Write("", TObject::kOverwrite);
		sinFitCalib->Write("", TObject::kOverwrite);

	}

	// </SALVATAGGIO>

	// </OUTPUT>

	delete sinFitCalib, sinEven, sinOdd;
	delete line3, line4;
	delete gr1, grEven, grOdd, grCalib;
	delete [] iSamp;
	delete [] ch;
	delete [] res;

	tf1.Close();

	return 0;

}


		//residui canale (double) - fit
		// c1->cd(6);
		// gPad->SetGrid();
		// gr3.SetTitle("Residui canale-fit in funzione del tempo");
		// gr3.GetXaxis()->SetTitle("t[s]");
		// gr3.GetYaxis()->SetTitle("#Delta[ch]");
		// gr3.DrawClone("PAL");
