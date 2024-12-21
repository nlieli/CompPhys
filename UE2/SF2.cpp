#include <iostream>
#include <fstream> // creates input file stream
#include <sstream>
#include <string>
#include <vector>
#include "SF2.h"

void readMatrix(const std::string& fileName, const delimiter& delim)
{
	std::ifstream fs;
	std::stringstream ss;
	std::string row;
	std::string value;
	int columns = 0;
	bool repeatedSpaces = false;

	fs.open(fileName, std::ios::in);
	if (!fs.is_open()) { std::cerr << "[ERROR] File stream did not open successfully!" << std::endl;  return; }

	ss << fs.rdbuf();
	std::getline(ss, row);

	if (delim == space)
	{
		for (size_t i = 0; i < row.size(); ++i)
		{
			if (row[i] == delim && !repeatedSpaces)
			{
				columns += 1;
				repeatedSpaces = true;
			}
			else
			{
				repeatedSpaces = false;
			}
		}
	} 
	else
	{
		for (size_t i = 0; i < row.size(); ++i)
		{
			if (row[i] == delim)
				columns += 1;
		}
	}

	std::vector<std::vector<double>> Matrix(columns); // use columns as rows because iterating over columns is faster and they typically contain more entries
	
	while (std::getline(ss, value, static_cast<char>(delim)))
	{
		size_t i = 0;
		
		Matrix[i % columns].push_back(2);
		i += 1;
	}
	
	
}
