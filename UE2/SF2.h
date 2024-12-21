#pragma once
#include <fstream>

enum decimalStyle
{
	comma = ',',
	dot = '.'
};

void readMatrix(const std::string& fileName, const decimalStyle& style = dot);

template <typename T> 
void printVector(std::vector<T>& vector)
{
	std::cout << '[';
	for (size_t i = 0; i < vector.size(); ++i)
		std::cout << vector[i] << ' ';
	std::cout << ']' << std::endl;
}
