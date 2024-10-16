#include <vector>
#include "VectorMath.h"
#include <iostream>
#include <cmath>
#include "SF.h"

double dx = 1e-4;

// functions labeled "polyXXX" are functions that return the coefficients of the polynomial in the order: [1, x, x^2, ...] 

Timer::Timer()
{
	start = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration<float> (0.0f);
	end = start;
}

Timer::~Timer()
{
	end = std::chrono::high_resolution_clock::now();
	duration = end - start;

	float ms = duration.count() * 1000.0f;
	std::cout << "Timer took " << ms << "ms " << std::endl;
}

unsigned long long int factorial(unsigned int x) // replace with gamma function ?? 
{
	unsigned long long int y = 1;

	for (int i = x; i > 1; --i)
		y *= i;

	return y;
}

double factorialCrazy(double x) // HOW ACCURATE IS THIS THING ?!??!?!??!
{
	double factorial = 1;
	if (x == 0)
		factorial = 1.0;
	else
		factorial = x * factorialCrazy(x - 1.0);
	return factorial;
}

unsigned long long int nCr(unsigned long long int n, unsigned long long int r) 
{
	return factorial(n) / (factorial(r) * factorial(n - r));
}

std::vector<double> polyLegendre(int n) 
{
	std::vector<double> Pn(n + 1);
	for (int k = 0; k <= n / 2; ++k)
		Pn[n - 2 * k] = std::pow(-1, k) * factorialCrazy(2 * n - 2 * k) / (factorialCrazy(n - k) * factorialCrazy(n - 2 * k) * factorialCrazy(k) * std::pow(2, n));
		//Pn[n - 2 * k] = std::pow(-1, k) * nCr(n, k) * nCr(2 * n - 2 * k, n) * std::pow(2, -n);
	return Pn;
}

std::vector<double> polyLegendreRec(int n) // recursive function
{
	std::vector<double> Pn(n + 1); // n-th legendre Polynomial	
	std::vector<double> Pn_(n + 1); // (n - 1)-th legendre Polynomial 
	std::vector<double> Pnx(n + 1); // (n + 1)-the legendre Polynomial 

	Pn_[0] = 1;
	Pn[0] = 0;
	Pn[1] = 1;

	if (n == 0) 
		return Pn_;

	if (n == 1)
		return Pn;

	for (int i = 1; i < n; ++i)
	{
		// Pnx = vectorDiv(vectorSub(vectorMultiply(2 * i + 1, vectorShiftRight(Pn)), vectorMultiply(Pn_, i)), i + 1);
		Pnx = (((2 * i + 1) * Pn >> 1) - (i * Pn_)) / (i + 1);
		Pn_ = Pn;
		Pn = Pnx;
	}

	return Pnx;
}

double polyNewtonRaphson(std::vector<double>& function, int N, double x0) // N = number of iterations, x0 = first guess
{
	for (int i = 0; i < N; ++i)
	{
		if (x0 == 0)
			return x0;
		x0 = x0 - (polyEval(function, x0) / polyEval(polyDiff(function), x0));
	}
		
	return x0;
}

std::vector<double> polyBracketing(std::vector<double>& function, double a, double b) // breaks down for higher orders - fix????
{
	std::vector<double> guess;
	for (double i = 0; (a + i) < b; i += dx)
	{
		if (polyBrackSearch(function, a + i, a + i + dx))
			guess.push_back(a + i + dx / 2);
	}

	return guess;
}

bool polyBrackSearch(std::vector<double>& function, double a, double b)
{
	double y1 = polyEval(function, a);
	double y2 = polyEval(function, b);
	double y = y1 * y2;
	if (y < 0)
		return true;

	return false;
}

std::vector<double> polyFindRoots(std::vector<double>& function, int N, double a, double b)
{
	std::vector<double> guess = polyBracketing(function, a, b);
	size_t numberOfRoots = guess.size();
	std::vector<double> result(numberOfRoots);

	for (size_t i = 0; i < numberOfRoots; ++i)
		result[i] = polyNewtonRaphson(function, N, guess[i]);

	return result;
}

std::vector<double> gaussLegendreWeights(std::vector<double> xi, std::vector<double> legendrePolynomial)
{
	size_t order = xi.size();
	std::vector<double> weights(order);

	for (size_t i = 0; i < order; ++i)
	{
		double PnPrime = polyEval(polyDiff(legendrePolynomial), xi[i]);
		weights[i] = 2 / ((1 - xi[i] * xi[i]) * PnPrime * PnPrime);
	}

	return weights;
}







