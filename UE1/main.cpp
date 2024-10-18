#include <iostream>
#include <vector>
#include "VectorMath.h" // functions for operations applied to vectors
#include <cmath> // std::sqrt, pow and other math functions
#include <iomanip> // std::setprecision 
#include <string_view>
#include <direct.h>
#include "matplot/matplot.h" // maybe plotting 

#include "SF.h" // various supporting functions

static double function(double x)
{
	return 1 / std::sqrt(1 - std::pow( x, 4));
}

static double trapezoid(double(*integrand)(double x), int N, double a, double b) 
{
	Timer timer;
	double dx = (b - a) / N;
	double F = 0;

	for (int i = 1; i < N; ++i)
	{
		F += dx / 2 * (integrand(a) + integrand(a + dx));
		a += dx;
	}

	return F;
}

static double simpson(double(*integrand)(double x), int N, double a, double b)
{
	Timer timer;
	double dx = (b - a) / N;
	double F = 0;
	double S1 = 0;
	double S2 = 0;
	double a0 = a;
	
	for (int i = 1; i < N; ++i)
	{
		a += dx;
		if (i % 2)
			S1 += integrand(a); // sum over odd indices	
		else
			S2 += integrand(a); // sum over even indices
	}

	F = dx / 3 * (integrand(a0) + integrand(a) + 4 * S1 + 2 * S2); // Simpson's Rule (integrand(a) sollte integrand(b) = inf sein)
	
	return F;
}

static double gaussquad(double(*integrand)(double x), int N, double a, double b)
{
	Timer timer;
	double F = 0;
	double chi;
	std::vector<double> w(N);
	std::vector<double> x(N);
	std::vector<double> P(N + 1); // to account for x^0;

	P = polyLegendreRec(N); // returns vector with coefficients for legendre polynomial in order [1, x, x^2, ... , x^N]
	x = polyFindRoots(P, N, -1, 1); // boundaries -1, 1 because Legendre polynomials always have all their roots in between those boundaries
	w = gaussLegendreWeights(x, P);  // weights for Gauss Legendre quadrature (int_a^b f(x) dx = sum_{i=1}^n wi f(xi));

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
	double b = 1; // upper bound for integration 
	int N = 0; // number of integration range subdivisions
	double Ft = trapezoid(function, 40, 0, 1);
	double Fs = simpson(function, 40, 0, 1);
	double Fg = gaussquad(function, 36, 0, 1); // for N > 40 polyBracketing function gives wrong values
	
	std::cout << std::setprecision(15) << "Trap = " << Ft << " Simp = " << Fs << " gauss = " << Fg << std::endl;

	//const char* gnuplotPath = "PATH=C:\\CPProjects\\CompPhys\\Dependencies\\gnuplot\\bin";
	const char* gnuplotPath = "PATH=..\\..\\Dependencies\\gnuplot\\bin";
	if (_putenv(gnuplotPath) == 0)
		std::cout << "Environment variable set successfully." << std::endl;
	else
		std::cerr << "Failed to set environment variable." << std::endl;



	std::vector<double> y = { 1, 2, 3, 4 };
	std::vector<double> x = y;

	matplot::plot(x, y);
	matplot::show();


}
