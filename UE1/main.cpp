// general purpose
#include <iostream>
#include <iomanip> // std::setprecision 
#include <direct.h>
#include <type_traits>

// math related functions and libraries
#include <vector>
#include "VectorMath.h" // functions for operations applied to vectors
#include <cmath> // std::sqrt, pow and other math functions

// Graphing and plotting libraries
#include "matplot/matplot.h"
#include "plot.h"

// various supporting functions
#include "SF.h" 
#include "recursive.h"

// Timer in debug mode if needed
#ifdef _DEBUG
#define TIMER Timer timer
#else 
#define TIMER
#endif


static double function(double x)
{
	return 1 / std::sqrt(1 - std::pow( x, 4));
}

static double T_integrand(double(*potential)(double x), double a, double x)
{
	return 1 / std::sqrt(potential(a) - potential(x));
}

static double trapezoid(double(*integrand)(double x), int N, double a, double b) 
{
	double dx = (b - a) / N;
	double F = 0;

	for (int i = 1; i <= N; ++i)
	{
		F += dx / 2 * (integrand(a) + integrand(a + dx));
		a += dx;
	}

	return F;
}


static double simpson(double(*integrand)(double x), int n, double a, double b)
{
	double dx = (b - a) / n;
	double f = 0;
	double S1 = 0;
	double S2 = 0;
	double a0 = a;
	
	for (int i = 1; i < n; ++i)
	{
		a += dx;
		if (i % 2)
			S1 += integrand(a); // sum over odd indices	
		else
			S2 += integrand(a); // sum over even indices
	}

	f = dx / 3 * (integrand(a0) + integrand(a) + 4 * S1 + 2 * S2); // simpson's rule (integrand(a) sollte integrand(b) = inf sein)
	
	return f;
}

static double legendreGauss(double(*integrand)(double x), int N, double a, double b)
{
	double F = 0;
	double chi;
	std::vector<double> w(N);	
	std::vector<double> x(N);
	std::vector<double> P(N + 1); // to account for x^0;
	
	x = legendreFindRoots(N, 100, -1, 1);
	w = legendreWeights(x, N);

	for (int i = 0; i < N; ++i)
	{
		chi = ((b - a) / 2) * x[i] + ((a + b) / 2);
		F += w[i] * integrand(chi); // integrand boundaries get shifted from a,b to [-1, 1] to fit Gauss Quadrature rule
	}
	F = (b - a) / 2 * F;

	return F;
}

int main() 
{ // compare to which value?
	double a = 0; // lower bound for integration
	double b = 0.999; // upper bound for integration 
	int N = 100; // number of integration range subdivisions
	double trueValue = 1.311028777146120;
	
	std::vector<std::vector<double>> results(3, std::vector<double>(N));
	
	for (int i = 1; i < N - 1; ++i)
	{
		results[0][i - 1] = trapezoid(function, i, a, b);
		results[2][i - 1] = legendreGauss(function, i, a, b); // for N > 40 polyBracketing function gives wrong values
	}	
	
	int j = 0;

	for (int i = 1; i < 2*N; i += 2)
	{
		results[1][j] = simpson(function, i, a, b);
		j += 1;
	}

	printVector(results[0]);
	printVector(results[1]);
	printVector(results[2]);

	T_integrand(cosh, 2, 3);
	T_integrand(exp, abs(a), abs(3));
	T_integrand([] (double x) {return -cos(x); }, a, 3);


	plotResults(results); 

	//int N = 350;
	//double f = legendreGauss(function, N, a, b);
	//double f = legendrePrimeEval(3, -0.774596);
	//std::cout << std::setprecision(15) << f << std::endl;

	//std::vector<double> x = matplot::linspace(-1, 1, 10000);
	//std::vector<double> y(x.size());
	//for (size_t i = 0; i < x.size(); ++i)
		//y[i] = legendreEval(1000, x[i]);

	//matplot::plot(x, y);
	//matplot::show();


}
