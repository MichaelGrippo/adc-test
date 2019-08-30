
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <stdlib.h>

using namespace std;

const bool debug = 0;
const int isInverted0 = 0;
const int isInverted1 = 0;
const int isInverted2 = 0;
const int isInverted3 = 0;

//function which takes the first 64 bits and looks for
//a shift in the data streaming.
//the shift will be used to skip N bits and realign all the channels

int LookingForShift(vector <vector <string> > data, int adc0column, unsigned long long pattern0, unsigned long long mask0) {

	bool found = 0;
	int shift0 = 0;

	string temp64bit = data[1][adc0column]+data[2][adc0column]+data[3][adc0column]+data[4][adc0column]+data[5][adc0column]+data[6][adc0column]+data[7][adc0column]+data[8][adc0column];
	if (debug) cout << "temp64bit: "<< temp64bit << endl;
	unsigned long long temp = stoull( temp64bit, nullptr, 2 );
	if (debug) cout << "temp: " << temp << endl;

	for ( int i = 44; i > -1; i-- ) {

		if (debug) cout << temp64bit << " hex: " << hex << temp << " i " << i << " i am doing the and with " << ( mask0 << (i) ) << " results "<< (  temp & ( mask0 << (i) )  ) << " looking for " << ( pattern0 << (i) ) << endl;

		if (   (  temp & ( mask0 << (i) )  ) == pattern0 << (i)   ) {

			if ( shift0 == 0 ) shift0 = 45-i;
			if ( debug ) cout << hex << "Pattern found: " << temp << " line " << " " << temp64bit << " shift " << shift0 << endl;
			found = 1;
			break;

		}

	} // I found the shift

	if ( found == 0 ) {

		cout << "pattern of channel not found" << endl;

		return -1;

	} else return shift0;

}

//Function to extract data from the file.
//it uses the shift calculated previously.
//It read 64 bits every time, but then it uses only 32 bits, skipping what it has already read and what has to be shifted.

void ExtractData( vector <vector <string> > data, int adc0column, unsigned long long pattern0, unsigned long long mask, int shift0,int nsamples, vector<int> &vect1,vector<int> &vect2) {

	unsigned long long maskAdc = 0x0FFF0;
	unsigned long long maskAdc1 = 0xFFF;

	for ( int irow = 1; irow < nsamples - 4; irow += 4 ) {
		string temp64bit = data[irow][adc0column] + data[irow+1][adc0column] + data[irow+2][adc0column] + data[irow+3][adc0column] + data[irow+4][adc0column] + data[irow+5][adc0column]  + data[irow+6][adc0column]  + data[irow+7][adc0column];
		if (debug) cout << "now I have read " << temp64bit << endl;
		unsigned long long temp = stoull(temp64bit, nullptr, 2);
		if (debug) cout << "temp "<< dec << temp << endl;
		temp = temp << (shift0-1);
		if (debug) cout << "after removing most significant bits " << temp << endl;
		temp = ( temp >> (shift0) );
		temp = ( temp >> (32-shift0) );
		if (debug) cout << "final value "<< temp << endl;
		unsigned long long adc1  = 0;
		unsigned long long adc   = 0;
		if (   (  temp & ( mask << (12) )  ) == pattern0 << (12)   ) {
			adc = (  temp & ( maskAdc<<(12) )  );
			adc = adc >> (12);
			adc = adc >> 4;
			adc1 = (temp & maskAdc1);
			if (debug) cout << dec << "row " << irow << hex << " Pattern found: " << temp << " " << temp64bit << " adc " << hex << adc << " adc1 "<< adc1 << " shift " << shift0 << endl;
		}
		//if (irow>10) break;

		vect1.push_back(adc1); //i+1
		vect2.push_back(adc); // i+3
	}

}

//Writing samples to txt file

int testcsv( string file = "Files/11_07_board2_adctm_sin500k_clk162.csv", string whichADC = "L") {

	vector <vector <string> > data;
	ifstream infile( file );
	vector <string> record;

	int nsamples = 0;
	int ncolumns = 0;
	bool isfirstline = 1;

	while (infile) {

		string s;
		if ( !getline( infile, s ) ) break;

		istringstream ss(s);

		while (ss) {

			string s;
			if ( !getline( ss, s, ',' ) ) break;
			record.push_back(s);

		}

		data.push_back( record );
		nsamples++;
		if ( isfirstline == 1 ) ncolumns = record.size();
		record.clear();
		isfirstline = 0;

	}

	if ( !infile.eof() ) {
		cerr << "Fooey! I did not find any file!\n";
		return -1;
	}


	cout << "nsamples found: " << nsamples << endl;
	cout << "ncolumns found: " << ncolumns << endl;

	string nametofind0("adc0"), nametofind1("adc1"),nametofind2("adc2"),nametofind3("adc3");

	int adc0column;
	int adc1column;
	int adc2column;
	int adc3column;

	//loop on the first line to check who is who
	for( int icol=0; icol < ncolumns; icol++ ) {

		if ( data[0][icol].find(nametofind0) != std::string::npos ) adc0column = icol;
		if ( data[0][icol].find(nametofind1) != std::string::npos ) adc1column = icol;
		if ( data[0][icol].find(nametofind2) != std::string::npos ) adc2column = icol;
		if ( data[0][icol].find(nametofind3) != std::string::npos ) adc3column = icol;

	}

	cout << "adc0 column: " << adc0column << endl;
	cout << "adc1 column: " << adc1column << endl;
	cout << "adc2 column: " << adc2column << endl;
	cout << "adc3 column: " << adc3column << endl;


	//extrapolate counts for each ADC:

	unsigned long long mask = 0xF000F;

	unsigned long long pattern0(0x30009), pattern1(0x6000c), pattern2(0xc0006), pattern3(0x90003);
	unsigned long long shift0(-1), shift1(-1), shift2(-1), shift3(-1);

	//I use the first line to check the shift:
	cout << "*****Looking for shift of channels" << endl;

	shift0 = LookingForShift(data,adc0column,pattern0,mask);
	shift1 = LookingForShift(data,adc1column,pattern1,mask);
	shift2 = LookingForShift(data,adc2column,pattern2,mask);
	shift3 = LookingForShift(data,adc3column,pattern3,mask);

	if ( shift0 == -1 ) cout << "PATTERN OF CHANNEL 0 NOT FOUND" << endl;
	if ( shift1 == -1 ) cout << "PATTERN OF CHANNEL 1 NOT FOUND" << endl;
	if ( shift2 == -1 ) cout << "PATTERN OF CHANNEL 2 NOT FOUND" << endl;
	if ( shift3 == -1 ) cout << "PATTERN OF CHANNEL 3 NOT FOUND" << endl;

	cout << "shift 0: " << shift0 << endl;
	cout << "shift 1: " << shift1 << endl;
	cout << "shift 2: " << shift2 << endl;
	cout << "shift 3: " << shift3 << endl;

	unsigned long long adc = 0;
	unsigned long long adc1 = 0;

	vector<int> Dout0_0;
	vector<int> Dout0_1;
	vector<int> Dout1_0;
	vector<int> Dout1_1;
	vector<int> Dout2_0;
	vector<int> Dout2_1;
	vector<int> Dout3_0;
	vector<int> Dout3_1;

	cout << "Restart after the shift: " << endl << endl;

	string temp64bit;

	//Extract data
	cout << "Extracting data from channel 0..." << endl;
	if ( shift0 != -1 ) ExtractData(data,adc0column,pattern0,mask,shift0,nsamples, Dout0_0, Dout0_1);
	cout << "Extracting data from channel 1..." << endl;
	if ( shift1 != -1 ) ExtractData(data,adc1column,pattern1,mask,shift1,nsamples, Dout1_0, Dout1_1);
	cout << "Extracting data from channel 2..." << endl;
	if ( shift2 != -1 ) ExtractData(data,adc2column,pattern2,mask,shift2,nsamples, Dout2_0, Dout2_1);
	cout << "Extracting data from channel 3..." << endl;
	if ( shift3 != -1 ) ExtractData(data,adc3column,pattern3,mask,shift3,nsamples, Dout3_0, Dout3_1);

	TH1D * ADC_H = new TH1D( "adcH", "adcH", 32700, -0.5, 32699.5 );
	TH1D * ADC_L = new TH1D( "adcL", "adcL", 32700, -0.5, 32699.5 );
	TH1D * ADC_L_Odd = new TH1D( "adcL_Odd", "adcL_Odd", 32700, -0.5, 32699.5 );
	TH1D * ADC_L_Even = new TH1D( "adcL_Even", "adcL_Even", 32700, -0.5, 32699.5 );

	TH1D * sampledistance = new TH1D( "sampledistance", "sampledistance", 400, -200, 200);

	int bin = 0;

	for ( int i = 0; i < Dout0_0.size(); i++ ) {
		ADC_H->SetBinContent( bin, Dout1_0[i] );
		ADC_H->SetBinContent( bin+1, Dout0_0[i] );
		ADC_H->SetBinContent( bin+2, Dout1_1[i] );
		ADC_H->SetBinContent( bin+3, Dout0_1[i] );
		bin += 4;
	}


	bin = 0;
	ofstream outfile;
	string outfilename = file + ".txt";
	cout << "ADC counts written in file " << outfilename << " for debugging" << endl;
	outfile.open(outfilename);

	int previoussample = 0;

	for( int i = 0; i < Dout2_0.size(); i++ ) {   //problema: usa solo ADC_L!

		ADC_L->SetBinContent( bin, Dout3_0[i] );
		ADC_L->SetBinContent( bin+1, Dout2_0[i] );
		ADC_L->SetBinContent( bin+2, Dout3_1[i] );
		ADC_L->SetBinContent( bin+3, Dout2_1[i] );

		ADC_L_Odd->SetBinContent( bin, Dout3_0[i] );
		ADC_L_Odd->SetBinContent( bin+2, Dout3_1[i] );

		ADC_L_Even->SetBinContent( bin+1, Dout2_0[i] );
		ADC_L_Even->SetBinContent( bin+3, Dout2_1[i] );

		//ADC LOW
		std::bitset<12> x1l( Dout3_0[i] );
		std::bitset<12> x2l( Dout2_0[i] );
		std::bitset<12> x3l( Dout3_1[i] );
		std::bitset<12> x4l( Dout2_1[i] );

		//ADC HIGH
		std::bitset<12> x1h( Dout1_0[i] );
		std::bitset<12> x2h( Dout0_0[i] );
		std::bitset<12> x3h( Dout1_1[i] );
		std::bitset<12> x4h( Dout0_1[i] );

		//SALVATAGGIO
		if ( whichADC != "H" && whichADC != "h" && whichADC != "L" && whichADC != "l" ){

			if (debug) {

				cout << "Warning: incorrect value for whichADC variable" << endl;
				cout << "\twhichADC = \"" << whichADC.c_str() << "\"" << endl;
				cout << "Defaulting to L" << endl;

			}

			whichADC = "L";

		}

		if ( whichADC == "L" || whichADC == "l" ) {

			outfile << x1l << endl;
			outfile << x2l << endl;
			outfile << x3l << endl;
			outfile << x4l << endl;

		} else if ( whichADC == "H" || whichADC == "h" ){

			outfile << x1h << endl;
			outfile << x2h << endl;
			outfile << x3h << endl;
			outfile << x4h << endl;

		}

		bin += 4;

		sampledistance->Fill( Dout3_0[i] - previoussample );
		sampledistance->Fill( Dout2_0[i] - Dout3_0[i] );
		sampledistance->Fill( Dout3_1[i] - Dout2_0[i] );
		sampledistance->Fill( Dout2_1[i] - Dout3_1[i] );
		previoussample = Dout2_1[i];
	}

	outfile.close();

	if (debug) {

		TCanvas *c0 = new TCanvas("All","All", 1024, 600);
		ADC_H->GetXaxis()->SetRangeUser(100, 400);
		ADC_H->Draw();

		TCanvas *c1 = new TCanvas("AllLow", "All Low", 1024, 600);
		ADC_L->SetMarkerStyle(20);
		ADC_L->GetXaxis()->SetRangeUser(0, 1000);
		ADC_L->Draw("p");

		TCanvas *c2 = new TCanvas("oddLow", "odd Low", 1024, 600);
		ADC_L_Odd->SetMarkerStyle(20);
		ADC_L_Odd->SetMarkerColor(kRed);
		ADC_L_Odd->GetXaxis()->SetRangeUser(100, 400);
		ADC_L_Odd->Draw("p");

		TCanvas *c3 = new TCanvas("evenLow", "even Low", 1024, 600);
		ADC_L_Even->SetMarkerStyle(20);
		ADC_L_Even->GetXaxis()->SetRangeUser(100, 400);
		ADC_L_Even->Draw("p");

		TCanvas *c4 = new TCanvas("mixLow", "mix Low", 1024, 600);
		ADC_L_Odd->Draw("p");
		ADC_L_Odd->GetYaxis()->SetRangeUser(0, 4200);
		ADC_L_Even->Draw("same p");

		TLegend *leg = new TLegend(0.1, 0.7, 0.3, 0.9);
		leg->AddEntry( ADC_L_Odd, "Odd", "p");
		leg->AddEntry( ADC_L_Even, "Even", "p");
		leg->Draw();

	}

	if ( shift0 != -1 ) {
		cout << "nMin Element Channel 0_0 = " << dec << *min_element( Dout0_0.begin(), Dout0_0.end() ) << endl;
		cout << "nMin Element Channel 0_1 = " << dec << *min_element( Dout0_1.begin(), Dout0_1.end() ) << endl << endl;
		cout << "nMax Element Channel 0_0 = " << dec << *max_element( Dout0_0.begin(), Dout0_0.end() ) << endl;
		cout << "nMax Element Channel 0_1 = " << dec << *max_element( Dout0_1.begin(), Dout0_1.end() ) << endl << endl;
	}

	if ( shift1 != -1 ) {    //era shift0!=-1 ---> bug!
		cout << "nMin Element Channel 1_0 = " << dec << *min_element( Dout1_0.begin(), Dout1_0.end() ) << endl;
		cout << "nMin Element Channel 1_1 = " << dec << *min_element( Dout1_1.begin(), Dout1_1.end() ) << endl << endl;
		cout << "nMax Element Channel 1_0 = " << dec << *max_element( Dout1_0.begin(), Dout1_0.end() ) << endl;
		cout << "nMax Element Channel 1_1 = " << dec << *max_element( Dout1_1.begin(), Dout1_1.end() ) << endl << endl;
	}

	if ( shift2 != -1 ) {
		cout << "nMin Element Channel 2_0 = " << dec << *min_element( Dout2_0.begin(), Dout2_0.end() ) << endl;
		cout << "nMin Element Channel 2_1 = " << dec << *min_element( Dout2_1.begin(), Dout2_1.end() ) << endl << endl;
		cout << "nMax Element Channel 2_0 = " << dec << *max_element( Dout2_0.begin(), Dout2_0.end() ) << endl;
		cout << "nMax Element Channel 2_1 = " << dec << *max_element( Dout2_1.begin(), Dout2_1.end() ) << endl << endl;
	}

	if ( shift3 != -1 ) {
		cout << "nMin Element Channel 3_0 = " << dec << *min_element( Dout3_0.begin(), Dout3_0.end() ) << endl;
		cout << "nMin Element Channel 3_1 = " << dec << *min_element( Dout3_1.begin(), Dout3_1.end() ) << endl << endl;
		cout << "nMax Element Channel 3_0 = " << dec << *max_element( Dout3_0.begin(), Dout3_0.end() ) << endl;
		cout << "nMax Element Channel 3_1 = " << dec << *max_element( Dout3_1.begin(), Dout3_1.end() ) << endl;
	}

	/*
	TH1 *hm =0;
	TVirtualFFT::SetTransform(0);
	hm = ADC_L->FFT(hm, "MAG");
	hm->SetTitle("Magnitude of the 1st transform");
	hm->Draw();
	*/

	return 0;

}
