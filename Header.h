#ifndef HEADER_H
#define HEADER_H

#include <string>
#include <unordered_map>

class Header : public TObject{

	public:

		Header();
		int AddInfo( std::string key, double val );
		int SetHeader( std::unordered_map<std::string, double> );
		std::unordered_map<std::string, double> GetHeader() const;

	private:

		std::unordered_map<std::string, double> header;

	ClassDef(Header, 1);

};

#endif