#pragma once
#include <vector>
#include "VectorMath.h"
#include <iostream>
#include <chrono>

// functions labeled "polyXXX" are functions that return the coefficients of the polynomial in the order: [1, x, x^2, ...] 

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

double polyNewtonRaphson(std::vector<double>& function, int N, double x0); // N = number of iterations, x0 = first guess

std::vector<double> polyFindRoots(std::vector<double>& function, int N, double a, double b);

std::vector<double> polyLegendre(int n); // n = n-th order legendre polynomial, 

std::vector<double> polyLegendreRec(int n); // n-th order legendre polynomial calculated recursively

std::vector<double> gaussLegendreWeights(std::vector<double> xi, std::vector<double> legendrePolynomial);

std::vector<double> polyBracketing(std::vector<double>& function, double a, double b); // function = polynomial in vector form, a = lower bound search area, b = upper bound search area

bool polyBrackSearch(std::vector<double>& function, double a, double b);


// -- templated functions

template <typename T>
struct is_vector : std::false_type {}; // helper functions that determine if input is of type std::vector<T>

template <typename T>
struct is_vector<std::vector<T>> : std::true_type {};

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
std::vector<T> polyDiff(std::vector<T> polynomial)
{
	static_assert(std::is_arithmetic<T>::value, "Polynomials must contain values of arithmetic type");
	size_t order = polynomial.size();
	for (size_t i = 1; i < order; ++i)
		polynomial[i - 1] = i * polynomial[i];

	polynomial[order - 1] = 0;

	return polynomial;
}







