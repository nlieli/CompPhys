#include <iostream>
#include <fstream> // creates input file stream
#include <sstream>
#include <string>
#include <regex>
#include <vector>
#include "SF2.h"

void readMatrix(const std::string& fileName, const decimalStyle& style)
{
	std::ifstream fs;
	std::stringstream ss;
	std::string row;
	std::vector<std::vector<double>> Matrix;
	std::regex reg;
	size_t columns = 0;

	if (style == dot)
		reg = "(-?\\d*\\.?\\d+e?[-+]\\d*)"; // decimal point system, captures all numbers using this style 
	else
		reg = "(-?\\d*,?\\d+e?[-+]\\d*)"; // decimal comma system, captures all numbers using this style

	fs.open(fileName, std::ios::in);
	if (!fs.is_open()) { std::cerr << "[ERROR] File stream did not open successfully!" << std::endl;  return; }

	ss << fs.rdbuf();
	std::getline(ss, row);
	std::smatch matches;
	std::regex_search(row, matches, reg);
	columns = matches.size();
	Matrix.reserve(columns);

	for (int i = 0; std::getline(ss, row); ++i)
	{
		for (int j = 0; std::regex_search(row, matches, reg); ++j)
		{
			Matrix[i].push_back(stod(matches.str())); // out of bounds error somewhere here
			row = matches.suffix().str();
		}
	}
	

}
