#pragma once
#include <vector>
#include <iostream>

unsigned int factorial(unsigned int x);

unsigned int nCr(unsigned int n, unsigned int r);

std::vector<double> nthlegendre(int n);

template <typename T> // make output more readable
void printVector(std::vector<T>& vector)
{
	std::cout << '[';
	for (size_t i = 0; i < vector.size(); ++i)
		std::cout << vector[i] << ' ';
	std::cout << ']' << std::endl;
}

template <typename T>
std::vector<T>& polydiff(std::vector<T>& polynomial)
{
	size_t order = polynomial.size();
	for (size_t i = 1; i < order; ++i)
	{
		polynomial[i - 1] = i * polynomial[i];
	}
	polynomial[order - 1] = 0;

	return polynomial;
}
