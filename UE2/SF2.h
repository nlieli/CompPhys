#pragma once
#include <fstream>

enum delimiter
{
	tab = '\t',
	space = ' ',
	comma = ',',
	semicolon = ';',
	colon = ':'
};

void readMatrix(const std::string& fileName, const delimiter& delim = space);

template <typename T> 
void printVector(std::vector<T>& vector)
{
	std::cout << '[';
	for (size_t i = 0; i < vector.size(); ++i)
		std::cout << vector[i] << ' ';
	std::cout << ']' << std::endl;
}
