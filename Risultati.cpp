
#include "Risultati.h"

#include <bitset>
#include <string>
#include <unordered_map>
#include <iostream>

ClassImp(Risultati);

//costruttore di default
Risultati::Risultati(){

	std::unordered_map<std::string, double> histDataTemp = {};
	std::unordered_map<std::string, double> fitDataTemp = {};
	std::unordered_map<std::string, double> fftDataTemp = {};

	this->histData = histDataTemp;
	this->fitData = fitDataTemp;
	this->fftData = fftDataTemp;

	std::bitset<3> tempBitMask(0);
	this->bitMask = tempBitMask;

}

int Risultati::AddHistData( std::string key, double val ){

	//~ for ( auto x : this->histData ) if ()
	this->histData[key] = val;
	if ( this->bitMask[0] == 0 ) this->bitMask[0] = 1;
	return 0;

}
int Risultati::AddFitData( std::string key, double val ){

	this->fitData[key] = val;
	if ( this->bitMask[1] == 0 ) this->bitMask[1] = 1;
	return 0;

}
int Risultati::AddFftData( std::string key, double val ){

	this->fftData[key] = val;
	if ( this->bitMask[2] == 0 ) this->bitMask[2] = 1;
	return 0;

}

//getter
std::unordered_map<std::string, double> Risultati::GetHistData() const { return this->histData; }
std::unordered_map<std::string, double> Risultati::GetFitData() const { return this->fitData; }
std::unordered_map<std::string, double> Risultati::GetFftData() const { return this->fftData; }
std::bitset<3> Risultati::GetBitMask() const { return this->bitMask; };

//setter
int Risultati::SetHistData( std::unordered_map<std::string, double> histData ) {

	if ( histData.size() != 0  ){

		this->bitMask[0] = 1;
		this->histData = histData;
		return 0;

	} else {

		std::cout << "unordered_map vuota; nulla da aggiungere" << std::endl;
		return -1;

	}

}
int Risultati::SetFitData( std::unordered_map<std::string, double> fitData ) {

	if ( fitData.size() != 0  ){

		this->bitMask[1] = 1;
		this->fitData = fitData;
		return 0;

	} else {

		std::cout << "unordered_map vuota; nulla da aggiungere" << std::endl;
		return -1;

	}

}
int Risultati::SetFftData( std::unordered_map<std::string, double> fftData ) {

	if ( fftData.size() != 0  ){

		this->bitMask[2] = 1;
		this->fftData = fftData;
		return 0;

	} else {

		std::cout << "unordered_map vuota; nulla da aggiungere" << std::endl;
		return -1;

	}

}

