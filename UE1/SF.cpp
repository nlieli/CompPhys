#include <vector>
#include <iostream>
#include <cmath>

unsigned int factorial(unsigned int x) // should be generalized
{
	int y = 1;

	for (int i = x; i > 1; --i)
		y *= i;

	return y;
}

std::vector<double>& nthlegendre(size_t n)
{
	std::vector<double> Pn(n);
	size_t l;

	if (n % 2)
		l = (n - 1) / 2; // if n is odd
	else
		l = n / 2; // if n is even

	for (size_t k = 0; k < l; ++k)
		Pn[k] = std::pow(-1, k) * factorial((2 * n - 2 * k)) / (factorial(n - k) * factorial(n - 2 * k) * factorial(k) * std::pow(2, n)); // WRONG!!! 0 terms!

	return Pn;
}

template <typename T>
void printVector(std::vector<T>& vector)
{
	for (size_t i = 0; i < vector.size(); ++i)
		std::cout << vector[i];

	std::cout << std::endl;

}