
//~ #include "Header.h"

#if !defined (__CLING__) || defined (__ROOTCLING__)
#include <TH1D.h>
#include <TF1.h>
#include <TMath.h>
#include <math.h>
#include <string>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TCanvas.h>
#include <TPad.h>
#include <fstream>
#include <iostream>
#include <string>
#include <cstdlib>		//per usare atof
#include <vector>
#include <bitset>
#endif

void txtToROOT(

	std::string nomeFile = "Dati/dati.txt",
	bool save = 0,
	bool print = 1,
	bool binaryMode = 1,
	bool split = 0		//divide i campionamenti pari e dispari

	){

	//definizioni varie
	std::string nomeIst = "histSamp";
	std::string nomeIstEven = nomeIst + "Even";
	std::string nomeIstOdd = nomeIst + "Odd";
	std::string nomeTGraph = "grSamp";
	std::string nomeTGraphEven = nomeTGraph + "Even";
	std::string nomeTGraphOdd = nomeTGraph + "Odd";
	std::string nomeOutput = nomeFile;
	for (int i = 0; i < 4; i++) nomeOutput.pop_back();		//assumendo ".txt"
	nomeOutput.append(".root");

	//definizioni iniziali parametri ADC
	int nBit = 12;
	int nCh = (int)TMath::Power(2, nBit);
	int chMax = nCh - 1;

	//istogramma frequenze
	TH1D hist( nomeIst.c_str(), "Distribution of ADC codes", nCh, -0.5, chMax + 0.5 );
	hist.GetXaxis()->SetTitle("ADC code");
	hist.GetYaxis()->SetTitle("Frequency");
	hist.SetFillColor(kRed + 2);
	hist.SetLineColor(kRed + 2);

	//istogramma frequenze, campionamenti pari
	TH1D histEven( nomeIstEven.c_str(), "Distribution of ADC codes, even samples", nCh, -0.5, chMax + 0.5 );
	histEven.GetXaxis()->SetTitle("ADC code");
	histEven.GetYaxis()->SetTitle("Frequency");
	histEven.SetFillColor(kRed + 2);
	histEven.SetLineColor(kRed + 2);

	//istogramma frequenze, campionamenti dispari
	TH1D histOdd( nomeIstOdd.c_str(), "Distribution of ADC codes, odd samples", nCh, -0.5, chMax + 0.5 );
	histOdd.GetXaxis()->SetTitle("ADC code");
	histOdd.GetYaxis()->SetTitle("Frequency");
	histOdd.SetFillColor(kRed + 2);
	histOdd.SetLineColor(kRed + 2);

	//apro il file di testo
	std::ifstream fileDati;
	fileDati.open( nomeFile.c_str() );

	std::string riga = "";

	long i = 0;
	long iEven = 0;
	long iOdd = 0;
	int adc = 0;

	TGraph tGr;
	tGr.SetName( nomeTGraph.c_str() );
	tGr.SetTitle("ADC codes vs sample number");
	tGr.SetMarkerStyle(25);
	tGr.GetXaxis()->SetTitle("Sample number");
	tGr.GetYaxis()->SetTitle("ADC code");

	TGraph tGrEven;
	tGrEven.SetName( nomeTGraphEven.c_str() );
	tGrEven.SetTitle("ADC codes vs even sample number");
	tGrEven.SetMarkerStyle(25);
	tGrEven.GetXaxis()->SetTitle("Sample number");
	tGrEven.GetYaxis()->SetTitle("ADC code");

	TGraph tGrOdd;
	tGrOdd.SetName( nomeTGraphOdd.c_str() );
	tGrOdd.SetTitle("ADC codes vs odd sample number");
	tGrOdd.SetMarkerStyle(25);
	tGrOdd.GetXaxis()->SetTitle("Sample number");
	tGrOdd.GetYaxis()->SetTitle("ADC code");

	Header header;

	//ciclo sul file di testo
	if ( fileDati.is_open() ){

		while ( std::getline(fileDati, riga) ){

			int stringSize = riga.size();

			if ( riga.find("#") < stringSize ) {

				//estrazione key e value
				while ( riga.find(" ") < stringSize ) riga.erase( riga.find(" ") , 1 );	//rimozione degli spazi
				riga.erase( riga.find("#") , 1 );

				int pos = riga.find("=");
				std::string key = riga.substr(0, pos);
				std::string valString = riga.substr( pos + 1, stringSize - pos );
				double value = std::stod(valString);

				header.AddInfo( key, value );

				continue;

			}

			double adcRaw = -1;

			try{

				adcRaw = std::stod( riga );

			} catch ( const std::invalid_argument& invArg ) {

				std::cerr << adcRaw << " ---> " << "Argomento invalido a " << invArg.what() << " nella riga " << i + 1 << std::endl;
				continue;

			}

			if (!binaryMode) {

				adc = round( adcRaw );

			} else {

				adc = (int) std::strtoul( riga.c_str(), nullptr, 2 );

			}

			tGr.SetPoint( i, i, adc );
			hist.Fill( adc );

			if ( i % 2 == 0 ){

				tGrEven.SetPoint( iEven, iEven, adc );
				histEven.Fill( adc );

				iEven++;

			} else {

				tGrOdd.SetPoint( iOdd, iOdd, adc );
				histOdd.Fill( adc );

				iOdd++;

			}

			i++;	//<==== controllare questo

		}

	}

	const int nPts = i;
	const int nEven = iEven;
	const int nOdd = iOdd;
	fileDati.close();

	if (print) {

		//rette per campionamenti
		TF1	* line1 = new TF1( "line1", "0", 0, nPts - 1 );
		line1->SetLineColor(kOrange - 3);
		line1->SetLineWidth(1);
		line1->SetLineStyle(1);
		TF1	* line2 = new TF1( "line2", std::to_string(nCh).c_str(), 0, nPts - 1 );
		line2->SetLineColor(kOrange - 3);
		line2->SetLineWidth(1);
		line2->SetLineStyle(1);

		//grafica
		TCanvas * c1 = new TCanvas("c1");
		c1->Divide(2,1);
		c1->GetPad(1)->SetLeftMargin(0.15);			//serve per far sì che le label degli assi Y non
		c1->GetPad(2)->SetLeftMargin(0.15);			//vengano tagliate perché il margine è troppo piccolo

		c1->cd(1);
		gPad->SetGrid();
		hist.DrawCopy();

		c1->cd(2);
		gPad->SetGrid();
		tGr.SetMarkerStyle(20);
		tGr.SetMarkerSize(0.5);
		tGr.GetYaxis()->SetRangeUser(-100, nCh + 100);
		tGr.DrawClone("AP");
		line1->DrawCopy("same");
		line2->DrawCopy("same");

		delete line1;
		delete line2;

	}

	if (save){

		//salvataggio
		TFile tf1( nomeOutput.c_str(), "recreate");
		hist.Write();
		tGr.Write();

		if ( split ){

			tGrEven.Write();
			tGrOdd.Write();

			histEven.Write();
			histOdd.Write();

		}

		header.Write();
		tf1.Close();

	}

}
