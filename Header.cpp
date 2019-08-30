
#include "Header.h"

#include <string>
#include <unordered_map>
#include <iostream>

ClassImp(Header);

Header::Header(){

	std::unordered_map<std::string, double> tempMap = {};
	this->header = tempMap;

}

int Header::AddInfo( std::string key, double val ){

	this->header[key] = val;
	return 0;

}

std::unordered_map<std::string, double> Header::GetHeader() const {

	return this->header;

}

int Header::SetHeader( std::unordered_map<std::string, double> header ){

	if ( header.size() != 0  ){

		this->header = header;
		return 0;

	} else {

		std::cout << "unordered_map vuota; nulla da aggiungere" << std::endl;
		return -1;

	}

}