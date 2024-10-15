#pragma once
#include <vector>
#include <iostream>
#include <type_traits> // for type safety

// functions labeled "polyXXX" are functions that return the coefficients of the polynomial in the order: [1, x, x^2, ...] 

// --- general 

unsigned long long int factorial(unsigned int x);

unsigned long long int nCr(unsigned long long int n, unsigned long long int r);

// --- main SF.cpp function

double polyNewtonRaphson(std::vector<double>& function, int N, double x0); // N = number of iterations, x0 = first guess

std::vector<double> polyFindRoots(std::vector<double>& function, int N, double a, double b);

std::vector<double> polyLegendre(int n); // n = n-th order legendre polynomial, 

std::vector<double> gaussLegendreWeights(std::vector<double> xi, std::vector<double> legendrePolynomial);

std::vector<double> polyBracketing(std::vector<double>& function, double a, double b);

bool polyBrackSearch(std::vector<double>& function, double a, double b);


// -- templated functions

template <typename T> 
void printVector(std::vector<T>& vector)
{
	std::cout << '[';
	for (size_t i = 0; i < vector.size(); ++i)
		std::cout << vector[i] << ' ';
	std::cout << ']' << std::endl;
}

template <typename A, typename B>
double polyEval(std::vector<A> polynomial, B x)
{
	static_assert(std::is_arithmetic<A>::value && std::is_arithmetic<B>::value, "Polynomials must contain values of arithmetic type");
	double result = 0;
	size_t order = polynomial.size();
	for (size_t i = 0; i < order; ++i)
	{
		result += polynomial[i] * std::pow(x, i);
	}
	
	return result;
}

template <typename T>
std::vector<T> polyDiff(std::vector<T> polynomial)
{
	static_assert(std::is_arithmetic<T>::value, "Polynomials must contain values of arithmetic type");
	size_t order = polynomial.size();
	for (size_t i = 1; i < order; ++i)
		polynomial[i - 1] = i * polynomial[i];

	polynomial[order - 1] = 0;

	return polynomial;
}

template <typename T>
std::vector<T> vectorAdd(std::vector<T> v1, std::vector<T> v2) // Addition of two vectors
{
	static_assert(v1.size() == v2.size(), "Vectors must be the same length");
	size_t length = v1.size();
	for (size_t i = 0; i < length; ++i)
		v1[i] += v2[i];
	return v1;
}

template <typename T>
std::vector<T> operator+(std::vector<T> v1, std::vector<T> v2) // Overload for vector addition	
{
	return vectorAdd(v1, v2);
}

template <typename T>
std::vector<T> vectorSub(std::vector<T> v1, std::vector<T> v2) // Subtraction of two vectors follows common order 
{
	static_assert(v1.size() == v2.size(), "Vectors must be the same length");
	size_t length = v1.size();
	for (size_t i = 0; i < length; ++i)
		v1[i] -= v2[i];
	return v1;
}

template <typename T>
std::vector<T> operator-(std::vector<T> v1, std::vector<T> v2) // Overload for vector subtraction
{
	return vectorSub(v1, v2);
}

template <typename T>
std::vector<T> vectorMultiply(std::vector<T> vector, T scalar)
{
	size_t length = vector.size();
	for (size_t i = 0; i < length; ++i)
		vector[i] *= scalar;

	return vector;
}

template <typename T>
std::vector<T> vectorMultiply(T scalar, std::vector<T> vector)
{
	size_t length = vector.size();
	for (size_t i = 0; i < length; ++i)
		vector[i] *= scalar;

	return vector;
}

template <typename T>
std::vector<T> operator*(std::vector<T> vector, T scalar)
{
	return vectorMultiply(vector, scalar);
}

template <typename T>
std::vector<T> operator*(T scalar, std::vector<T> vector)
{
	return vectorMultiply(scalar, vector);
}

template <typename T>
std::vector<T> vectorDiv(std::vector<T> vector, T scalar)
{
	size_t length = vector.size();
	for (size_t i = 0; i < length; ++i)
		vector[i] /= scalar;
	
	return vector;
}

template <typename T>
std::vector<T> operator/(std::vector<T> vector, T scalar)
{
	return vectorDiv(vector, scalar);
}

template <typename T>
std::vector<T> vectorShiftRight(std::vector<T> vector)
{
	size_t length = vector.size();
	if (length == 0) return vector;

	for (size_t i = length; i > 0; --i)
		vector[i] = vector[i - 1];

	vector[0] = 0;
	return vector;
}

template <typename T>
std::vector<T> operator>>(std::vector<T> vector, int)
{
	return vectorShiftRight(vector);
}

template <typename T>
std::vector<T> vectorShiftLeft(std::vector<T> vector)
{
	size_t length = vector.size();
	for (size_t i = 1; i < length; ++i)
		vector[i - 1] = vector[i];

	vector[length - 1] = 0;
}

template <typename T>
std::vector<T> operator<<(std::vector<T> vector, int)
{
	return vectorShiftLeft(vector);
}





