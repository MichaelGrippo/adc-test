#ifndef RISULTATI_H
#define RISULTATI_H

#include <bitset>
#include <unordered_map>
#include <string>

//serve per salvare i dati nel TFile
class Risultati: public TObject{

	public:

		Risultati();

		int AddHistData( std::string key, double val );
		int AddFitData( std::string key, double val );
		int AddFftData( std::string key, double val );
		
		int SetHistData( std::unordered_map<std::string, double> );
		int SetFitData( std::unordered_map<std::string, double> );
		int SetFftData( std::unordered_map<std::string, double> );

		std::unordered_map<std::string, double> GetHistData() const;
		std::unordered_map<std::string, double> GetFitData() const;
		std::unordered_map<std::string, double> GetFftData() const;
		std::bitset<3> GetBitMask() const;

	private:

		std::unordered_map<std::string, double> histData;
		std::unordered_map<std::string, double> fitData;
		std::unordered_map<std::string, double> fftData;

		std::bitset<3> bitMask;

	ClassDef(Risultati, 1);

};

#endif