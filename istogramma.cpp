
#include "Risultati.h"
#include "Header.h"

#if !defined (__CLING__) || defined (__ROOTCLING__)
//		C++
#include <iostream>
#include <algorithm>
#include <unordered_map>
#include <string>
#include <cmath>
//		ROOT
#include <TH1D.h>
#include <TF1.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TStyle.h>
#include <TFile.h>
#endif

//	<DA FARSI>
//	-correggere cumulHist = *temp, cumulHistTh = *temp ed eventualmente altri casi
//	-fare controllo sulla buona riuscita del fit ===> fatto ma non soddisfacente
//	-controllare se bisogni usare m_T in caso di overflow ===> pare di no
//	</DA FARSI>

int istogramma(

	std::string nomeFile = "5e10-puro.root",
	bool salva = 0,
	bool grafica = 1,
	bool output = 1,
	bool useCalibratedSamples = 1,		//analizza i dati calibrati, se presenti
	int split = 0 		//0: tutti, 1: dispari, 2: pari

	){

	//definizioni iniziali parametri ADC
	const int nBit = 12;
	const double vmin = -0.6;
	const double vmax = +0.6;
	const int nCh = (int)TMath::Power(2, nBit);
	const int chMax = nCh - 1;
	const double LSB = (vmax - vmin)/nCh;

	//definizioni iniziali varie
	const double fitThreshold = 0.1;
	const std::string nomeIstDefault = "histSamp";
	std::string nomeIst = nomeIstDefault;
	if ( split == 1 && useCalibratedSamples == 0 ) nomeIst += "Odd";
	if ( split == 2 && useCalibratedSamples == 0 ) nomeIst += "Even";
	if ( useCalibratedSamples ) nomeIst += "Calib";
	std::string sampleType = "";
	if ( !useCalibratedSamples && !split ) sampleType += "All";
	if ( split == 1 ) sampleType += "Odd";
	if ( split == 2 ) sampleType += "Even";
	if ( useCalibratedSamples ) sampleType += "Calibrated";

	//~ const std::string nomeIst = "histSamples";
	const std::string ampKey = "amp[V]";
	const std::string offsetKey = "offset[V]";

	//retta
	TF1	* line = new TF1("line1", "0", 0, 4095);
	line->SetLineColor(kOrange - 3);
	line->SetLineWidth(1);
	line->SetLineStyle(1);

	//retta per fit
	TF1 * retta = new TF1("rettaTransizioni", "[0]*x + [1]", 1, chMax - 1);
	retta->SetParName(0, "m");
	retta->SetParName(1, "q");
	retta->SetLineColor(kOrange - 3);
	retta->SetLineWidth(1);
	retta->SetLineStyle(1);

	//retta per fit, valori ideali
	TF1 * rettaIdeale = new TF1("rettaTransizioniIdeali", "[0]*x + [1]", 1, chMax - 1);
	rettaIdeale->SetParName(0, "m");
	rettaIdeale->SetParName(1, "q");

	//	<DEFINIZIONI ISTOGRAMMI>

	//istogramma frequenze
	TH1D hist("histFreq", "Distribution of ADC codes", nCh, -0.5, chMax + 0.5);
	if (split == 1) hist.SetTitle("Distribution of ADC codes, odd samples");
	if (split == 2) hist.SetTitle("Distribution of ADC codes, even samples");
	hist.GetXaxis()->SetTitle("ADC code");
	hist.GetYaxis()->SetTitle("Counts");
	hist.SetFillColor(kRed + 2);
	hist.SetLineColor(kRed + 2);

	//istogramma frequenze ideale
	TH1D histADCth("histADCth", "Ideal distribution of ADC codes", nCh, -0.5, chMax + 0.5);
	histADCth.GetXaxis()->SetTitle("ADC code");
	histADCth.GetYaxis()->SetTitle("Counts");
	histADCth.SetLineColor(kOrange - 3);

	//istogramma cumulativo
	TH1D cumulHist("cumulHist", "Cumulative histogram", nCh, -0.5, chMax + 0.5);
	if (split == 1) cumulHist.SetTitle("Cumulative histogram, odd samples");
	if (split == 2) cumulHist.SetTitle("Cumulative histogram, even samples");
	cumulHist.GetXaxis()->SetTitle("ADC code");
	cumulHist.GetYaxis()->SetTitle("Cumulative counts");
	cumulHist.SetLineColor(kRed + 2);

	//istogramma cumulativo teorico
	TH1D cumulHistTh("cumulHistTh", "Cumulative histogram", nCh, -0.5, chMax + 0.5);
	cumulHistTh.GetXaxis()->SetTitle("ADC code");
	cumulHistTh.GetYaxis()->SetTitle("Cumulative counts");
	cumulHistTh.SetLineColorAlpha(kMagenta - 2, 0.5);

	//transizioni
	TH1D histTran("histTran", "Transition voltages", nCh, -0.5, chMax + 0.5);
	histTran.GetXaxis()->SetTitle("ADC code");
	histTran.GetYaxis()->SetTitle("Transition voltage[V]");
	histTran.SetLineColor(kRed + 2);

	//transizioni ideali
	TH1D histTranTh("histTranTh", "Ideal transition voltages", nCh, -0.5, chMax + 0.5);
	histTranTh.GetXaxis()->SetTitle("ADC code");
	histTranTh.GetYaxis()->SetTitle("Transition voltage[V]");
	histTranTh.SetLineColor(kOrange - 3);
	histTranTh.SetLineStyle(1);

	//transizioni corrette dopo il fit
	TH1D histTranFit("histTranFit", "Transition voltages after gain and offset correction", nCh, -0.5, chMax + 0.5);
	histTranFit.GetXaxis()->SetTitle("ADC count");
	histTranFit.GetYaxis()->SetTitle("V[V]");
	histTranFit.SetLineColor(kRed + 2);
	histTranFit.SetLineStyle(1);

	//DNL
	TH1D histDNL("histDNL", "Differential nonlinearity", nCh - 2, 0.5, chMax - 0.5);
	histDNL.GetXaxis()->SetTitle("ADC code");
	histDNL.GetYaxis()->SetTitle("DNL[LSB]");
	histDNL.GetYaxis()->SetRangeUser(-1, 1);
	histDNL.SetLineColor(kRed + 2);

	//DNL, no fit
	TH1D histDNL2("histDNL2", "Differential nonlinearity", nCh - 2, 0.5, chMax - 0.5);
	histDNL2.GetXaxis()->SetTitle("ADC count");
	histDNL2.GetYaxis()->SetTitle("DNL[LSB]");
	histDNL2.GetYaxis()->SetRangeUser(-1, 1);
	histDNL2.SetLineColor(kMagenta - 2);

	//INL
	TH1D histINL("histINL", "Integral nonlinearity", nCh, -0.5, chMax + 0.5);
	histINL.GetXaxis()->SetTitle("ADC code");
	histINL.GetYaxis()->SetTitle("INL[LSB]");
	histINL.GetYaxis()->SetRangeUser(-2, 2);
	histINL.SetLineColor(kRed + 2);

	//INL, no fit
	TH1D histINL2("histINL2", "Integral nonlinearity", nCh, -0.5, chMax + 0.5);
	histINL2.GetXaxis()->SetTitle("ADC code");
	histINL2.GetYaxis()->SetTitle("INL[LSB]");
	histINL2.GetYaxis()->SetRangeUser(-2, 2);
	histINL2.SetLineColor(kMagenta - 2);

	//	</DEFINIZIONI ISTOGRAMMI>

	//	<LETTURA DATI ADC>

	//apertura file ====> l'apertura del file ed il successivo controllo andrebbero messi all'inizio, però dà problemi
	TFile tf1( nomeFile.c_str(), "update");

	//controllo esistenza del file ====> sarebbe da mettere all'inizio
	if ( !tf1.IsOpen() ){

		std::cout << "File \"" << nomeFile.c_str() << "\" not found" << std::endl;
		std::cout << "Aborting execution" << std::endl;
		delete line, retta, rettaIdeale;
		return 1;

	}

	//ottiene l'istogramma dal file
	TH1D * histTemp = (TH1D *) tf1.Get( nomeIst.c_str() );

	if ( histTemp == nullptr && useCalibratedSamples ) {

		std::cout << "Histogram \"" << nomeIst.c_str() << "\" not found " << std::endl;
		std::cout << "Reverting to \"" << nomeIstDefault.c_str() << "\" " << std::endl;
		nomeIst = nomeIstDefault;
		histTemp = (TH1D *) tf1.Get( nomeIst.c_str() );

	}

	//controllo esistenza dell'istogramma
	if ( histTemp == nullptr ){

		std::cout << "Histogram \"" << nomeIst.c_str() << "\" not found" << std::endl;
		std::cout << "Aborting execution" << std::endl;
		delete line, retta, rettaIdeale, histTemp;
		tf1.Close();
		return 2;

	}

	//se esiste, lo passa in un altro TH1D
	hist = *histTemp;	//<==== aggiustare questo

	ULong_t m = hist.GetEntries();

	//lettura header
	Header * headerObj = (Header *) tf1.Get( "Header" );

	if ( headerObj == nullptr ){

		std::cout << "Header not found" << std::endl;
		std::cout << "Aborting execution" << std::endl;
		delete line, retta, rettaIdeale, histTemp;
		tf1.Close();
		return 3;

	}

	std::unordered_map < std::string, double > headerMap = headerObj->GetHeader();
	delete headerObj;

	double A = 0;
	double offset = 0;

	try {

		A = headerMap.at( ampKey );
		offset = headerMap.at( offsetKey );

	} catch ( const std::out_of_range& invArg ) {

		std::cout << "Key not found" << std::endl;
		std::cout << "Aborting execution" << std::endl;
		delete line, retta, rettaIdeale, histTemp;
		tf1.Close();
		return 4;

	}

	//~ tf1.Close();
	//~ delete histTemp;

	//	</LETTURA DATI ADC>

	//	<ANALISI>

	//ricavo alcuni parametri
	//~ ULong_t m_T = m - hist.GetBinContent(1) - hist.GetBinContent( chMax + 1 );
	//~ double A_exp = vmax / ( TMath::Sin( 0.5*TMath::Pi()*m_T/( m_T + hist.GetBinContent(1) + hist.GetBinContent(nCh) ) ) );

	//offset   <======

	//riempimento istogramma ADC teorico	===>	non considera un eventuale offset dell'onda. Non cambia nulla all'analisi ma correggilo
	for ( int k = 0; k < nCh + 1; k++ ){

		double va = ( k - 0.5 ) * LSB;
		double vb = ( k + 0.5 ) * LSB;

		double val = ( m/TMath::Pi() ) * (  TMath::ASin( (vb + vmin)/A ) - TMath::ASin( (va + vmin)/A )  );

		histADCth.SetBinContent( k + 1, val );

	}

	//calcolo istogramma cumulativo
	TH1D * temp = (TH1D *) hist.GetCumulative();
	cumulHist = *temp;		//correggere!
	temp = (TH1D *) histADCth.GetCumulative();
	cumulHistTh = *temp;	//correggere!

	//riempimento istogramma transizioni
	for ( int bin = 1; bin <= nCh; bin++ ){		//ciclo da ch 1 a ch 4095 (12 bit)

		// bin = 1 <-> ch = 0; bin = 2 <-> ch = 1; ...bin <-> ch - 1
		// bin + 1 <-> ch
		//~ double channel = bin - 1;
		histTran.SetBinContent(   bin + 1, offset - A * TMath::Cos(  ( TMath::Pi() / m ) * cumulHist.GetBinContent( bin )  )   );
		histTranTh.SetBinContent( bin + 1, ( bin - 0.5)*LSB + vmin );
		//~ histTranTh2.SetBinContent(   bin + 1, offset - A * TMath::Cos(  ( TMath::Pi() / m ) * cumulHistTh.GetBinContent( bin )  )   );

	}

	//		<FIT>

	//fit lineare alla caratteristica dell'ADC per ottenere errori di gain ed offset
	histTranTh.Fit("rettaTransizioniIdeali", "RQN0");
	//~ retta->SetParameter( 0, rettaIdeale->GetParameter(0) );
	//~ retta->SetParameter( 1, rettaIdeale->GetParameter(1) );
	histTran.Fit("rettaTransizioni", "RQN0");

	//conviene cambiare metodo: magari se m e q sono troppo diversi dai valori teorici
	if ( retta->GetChisquare() / retta->GetNDF() > fitThreshold ){

		std::cout << "WARNING: transition voltages fit did not converge" << std::endl;
		std::cout << "\tChi-square: " << retta->GetChisquare() << std::endl;
		std::cout << "\tDOF: " << retta->GetNDF() << std::endl;

	}


	double gainFit = retta->GetParameter("m");
	double offsetFit = retta->GetParameter("q");

	double gainFitErr = retta->GetParError(0);		//m
	double offsetFitErr = retta->GetParError(1);	//q

	double gainIdeale = rettaIdeale->GetParameter("m");
	double offsetIdeale = rettaIdeale->GetParameter("q");

	double G = gainFit / gainIdeale;
	double Vos = offsetFit - offsetIdeale;

	//		</FIT>

	//correzione delle tensioni di transizione con i gain e gli offset ricavati
	for (int bin = 2; bin <=  histTranTh.GetSize() - 2; bin ++){

		histTranFit.SetBinContent(  bin, ( histTran.GetBinContent(bin) - offsetFit )  / G + offsetIdeale  );	//controllare metodo: se si immette un'ampiezza sbagliata l'offset si sballa

	}

	//calcolo DNL
	for ( int ch = 1; ch < nCh - 1; ch ++ ){

		//con fit
		double val = ( histTranFit.GetBinContent( ch + 2 ) - histTranFit.GetBinContent( ch + 1 )  ) / LSB - 1;
		histDNL.SetBinContent( ch, val );

		//senza fit
		double val2 = ( histTran.GetBinContent( ch + 2 ) - histTran.GetBinContent( ch + 1 )  ) / LSB - 1;
		histDNL2.SetBinContent( ch, val2 );

	}

	//calcolo INL, metodo dello standard
	for( int bin = 2; bin <= histTranFit.GetSize() - 2 - 1; bin ++){	//-2 per bin di overflow, -1 per togliere l'ultimo bin

		//con fit
		double epsilon = histTranFit.GetBinContent( bin ) - histTranTh.GetBinContent( bin );
		histINL.SetBinContent( bin, epsilon / LSB );

		//senza fit
		double epsilon2 = histTran.GetBinContent( bin ) - histTranTh.GetBinContent( bin );
		histINL2.SetBinContent( bin, epsilon2 / LSB );

	}

	//missing codes
	std::vector<int> missingCodes;
	bool detectedMissingCodes = 0;

	for ( int bin = 1; bin <= histDNL.GetNbinsX(); bin++ ){

		if ( histDNL.GetBinContent(bin) < -0.9 ){

			missingCodes.push_back( bin - 1 );		//i bin partono da 1 e i canali da 0
			if (!detectedMissingCodes) detectedMissingCodes = 1;

		}

	}

	if (detectedMissingCodes) std::cout << "WARNING: missing codes detected" << std::endl;

	//	</ANALISI>

	//	<OUTPUT>

	if (output){

		std::cout << "=================== Results ===================" << std::endl;
		
		std::cout << "Dataset: " << nomeFile.c_str() << std::endl;
		std::cout << "Type: " << sampleType.c_str() << std::endl;

		std::cout << "Gain[V] = " << gainFit << std::endl;
		std::cout << "Fit error(gain)[V] = " << gainFitErr << std::endl;
		std::cout << "Ideal gain[V]= " << gainIdeale << std::endl;
		std::cout << "Offset[V] = " << offsetFit << std::endl;
		std::cout << "Fit error(offset)[V] = " << offsetFitErr << std::endl;
		std::cout << "Ideal offset[V] = " << offsetIdeale << std::endl;
		std::cout << "Gain / ideal gain = " << gainFit / gainIdeale << std::endl;
		std::cout << "Offset - ideal offset[V] = " << offsetFit - offsetIdeale << std::endl;

		std::cout << "DNL(max): " << histDNL.GetBinContent( histDNL.GetMaximumBin() ) << " LSB, code " << histDNL.GetMaximumBin() << std::endl;
		std::cout << "DNL(min): " << histDNL.GetBinContent( histDNL.GetMinimumBin() ) << " LSB, code " << histDNL.GetMinimumBin() << std::endl;
		std::cout << "INL(max): " << histINL.GetBinContent( histINL.GetMaximumBin() ) << " LSB, code " << histINL.GetMaximumBin() - 1 << std::endl;
		std::cout << "INL(min): " << histINL.GetBinContent( histINL.GetMinimumBin() ) << " LSB, code " << histINL.GetMinimumBin() - 1 << std::endl;

		std::cout << "Missing codes: " << missingCodes.size() << std::endl;
		//~ if (detectedMissingCodes) for ( int i = 0; i < missingCodes.size(); i++ ) std::cout << missingCodes.pop_back();
		
		std::cout << "===============================================" << std::endl;

	}

	if (grafica){

		//grafica
		TCanvas * c1 = new TCanvas("c1", "", 1250, 650);
		c1->Divide(3,2);
		for ( int i = 1; i < 7; i++ ) c1->GetPad( i )->SetLeftMargin(0.15);		//margini per non tagliare label Y
		//~ gStyle->SetOptStat(111);
		gStyle->SetStatX(0.75);
		gStyle->SetStatY(0.9);
		gStyle->SetOptStat(10);

		//istogramma conteggi
		c1->cd(1);
		gPad->SetGrid();
		hist.GetYaxis()->SetRangeUser( 0.95 * hist.GetMinimum(), 1.05 * hist.GetBinContent(2) );
		//~ hist.SetFillColorAlpha(kRed + 2, 0.05);
		hist.DrawCopy();
		histADCth.DrawCopy("same");

		//output ist. cumulativo
		c1->cd(2);
		gPad->SetGrid();
		cumulHist.SetFillColor(kWhite);
		cumulHistTh.SetFillColor(kWhite);
		cumulHistTh.SetLineColorAlpha(kOrange - 3, 0.2);
		cumulHist.GetYaxis()->SetRangeUser( 0, cumulHist.GetBinContent( cumulHist.GetMaximumBin() ) );
		cumulHist.SetTitle("Cumulative histogram");		//temp finché non trovo metodo migliore
		cumulHist.DrawCopy();
		cumulHistTh.DrawCopy("same");

		//output ist. transizioni
		c1->cd(3);
		gPad->SetGrid();
		histTranTh.SetTitle("Transition voltages");	//temp finché non trovo metodo migliore
		histTranTh.GetYaxis()->SetRangeUser(vmin, vmax);
		histTranTh.DrawCopy();
		//~ histTranTh2.DrawCopy("same");
		histTran.DrawCopy("same");
		retta->DrawCopy("same");
		gPad->Modified();

		//output correzione
		c1->cd(4);
		gPad->SetGrid();
		histTranTh.SetTitle("Transition voltages after gain and offset correction");
		histTranTh.DrawCopy();
		//~ histTranFit.GetYaxis()->SetRangeUser(vmin, vmax);
		histTranFit.DrawCopy("same");
		//~ histTran.DrawCopy("same");

		//output ist. DNL
		c1->cd(5);
		gPad->SetGrid();
		histDNL.DrawCopy();
		//~ histDNL2.DrawCopy("same");
		line->DrawCopy("same");

		//output ist INL
		c1->cd(6);
		gPad->SetGrid();
		histINL.SetTitle("Integral nonlinearity");
		histINL.DrawCopy();
		//~ histINL2.DrawCopy("same");
		line->DrawCopy("same");

		c1->Update();

		//~ c1->Print( ( nomeFile + "-output.pdf").c_str(), "pdf" );

	}

	//	</OUTPUT>

	//	<SALVATAGGIO>

	if (salva){

		double maxINL = 0;
		double maxDNL = 0;

		maxDNL = (   histDNL.GetBinContent(  histDNL.GetMaximumBin()  ) > histDNL.GetBinContent(  histDNL.GetMinimumBin()  )   ) ? histDNL.GetBinContent( histDNL.GetMaximumBin() ) : histDNL.GetBinContent( histDNL.GetMinimumBin() );
		maxINL = (   histINL.GetBinContent(  histINL.GetMaximumBin()  ) > histINL.GetBinContent(  histINL.GetMinimumBin()  )   ) ? histINL.GetBinContent( histINL.GetMaximumBin() ) : histINL.GetBinContent( histINL.GetMinimumBin() );

		Risultati * res = (Risultati *) tf1.Get("Risultati");
		if ( res == nullptr ) res = new Risultati();
		res->AddHistData( "DNL [LSB]", maxDNL );
		res->AddHistData( "INL [LSB]", maxINL );
		res->AddHistData( "Gain error [V]", gainFit );
		res->AddHistData( "Offset error [V]", offsetFit );
		res->AddHistData( "Gain error uncertainty [V]", gainFitErr );
		res->AddHistData( "Offset error uncertainty [V]", offsetFitErr );
		res->Write("", TObject::kOverwrite);
		delete res;

		hist.Write("", TObject::kOverwrite);
		cumulHist.Write("", TObject::kOverwrite);
		histTran.Write("", TObject::kOverwrite);
		histTranFit.Write("", TObject::kOverwrite);
		histDNL.Write("", TObject::kOverwrite);
		histINL.Write("", TObject::kOverwrite);

	}

	//	</SALVATAGGIO>

	//	<PULIZIA>

	delete line;
	delete retta, rettaIdeale;
	delete temp;

	tf1.Close();

	//	</PULIZIA>

	return 0;

}
