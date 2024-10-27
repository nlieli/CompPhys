#pragma once
#include <vector>
#include "VectorMath.h"
#include <iostream>
#include <chrono>

// functions labeled "polyXYZ" are functions that return the coefficients of the polynomial in the order: [1, x, x^2, ...] 

// --- general 

unsigned long long int factorial(unsigned int x);
double factorialCrazy(double x);

unsigned long long int nCr(unsigned long long int n, unsigned long long int r);

struct Timer
{
	std::chrono::time_point<std::chrono::steady_clock> start, end;
	std::chrono::duration<float> duration;

	Timer();

	~Timer();
};

// --- main SF.cpp function

double polyNewtonRaphson(std::vector<double>& function, int iterations, double x0); // x0 = first guess
double polyBisection(std::vector<double>& function, int iterations, std::array<double, 2> interval); // interval = rough interval in which the root is found 

std::vector<double> polyFindRoots(std::vector<double>& function, int iterations, double a, double b);

std::vector<double> polyLegendre(int n); // n = n-th order legendre polynomial, 
std::vector<double> polyLegendreRec(int n); // n-th order legendre polynomial calculated recursively

std::vector<double> gaussLegendreWeights(std::vector<double> xi, std::vector<double> legendrePolynomial);

std::vector<double> polyBracketing(std::vector<double>& function, double a, double b); // function = polynomial in vector form, a = lower bound search area, b = upper bound search area
std::vector<std::array<double, 2>> polyBracketingInterval(std::vector<double>& function, double a, double b); // polyBracketing function that outputs an interval instead of a guess

bool polyBrackSearch(const std::vector<double>& function, double a, double b);


// -- templated functions

template <typename T>
struct is_vector : std::false_type {}; // helper functions that determine if input is of type std::vector<T>

template <typename T>
struct is_vector<std::vector<T>> : std::true_type {};

template <typename T>
struct is_array: std::false_type {}; // helper functions that determine if input is of type std::array<T, size>

template <typename T, int size>
struct is_array<std::array<T, size>> : std::true_type {};

template <typename T> 
void printVector(std::vector<T>& vector)
{
	std::cout << '[';
	for (size_t i = 0; i < vector.size(); ++i)
		std::cout << vector[i] << ' ';
	std::cout << ']' << std::endl;
}

template <typename T, int size>
void printArray(std::array<T, size>& array)
{
	std::cout << '[';
	for (size_t i = 0; i < size; ++i)
		std::cout << array[i] << ' ';
	std::cout << ']' << std::endl;
}

template <typename T>
void print(T object) // function needs to be expanded to fit more types
{
	if constexpr (is_vector<T>::value) printVector(object);
	if constexpr (is_array<T>::value) printArray(object);
}

template <typename A, typename B>
double polyEval(std::vector<A> polynomial, B x) // evaluates polynomial of vector form at input position x
{
	static_assert(std::is_arithmetic<A>::value && std::is_arithmetic<B>::value, "Polynomials must contain values of arithmetic type");
	double result = 0;
	size_t order = polynomial.size();
	//for (size_t i = 0; i < order; ++i)
	//{
		//result += polynomial[i] * std::pow(x, i);
	//}

	for (size_t i = order; i > 0; --i)
	{
		result = result * x + polynomial[i - 1];
	}
	
	return result;
}

template <typename T>
std::vector<T> polyDiff(std::vector<T> polynomial) // ouputs the derivative of a polynomial function of vector form
{
	static_assert(std::is_arithmetic<T>::value, "Polynomials must contain values of arithmetic type");
	size_t order = polynomial.size();
	for (size_t i = 1; i < order; ++i)
		polynomial[i - 1] = i * polynomial[i];

	polynomial[order - 1] = 0;

	return polynomial;
}







