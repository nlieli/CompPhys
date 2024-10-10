#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <iomanip>

#include "SF.h"

static double function(double x)
{
	return 1 / std::sqrt(1 - std::pow( x, 4));
}

static double trapezoid(double(*integrand)(double x), int N, double a, double b)
{
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

static void gaussquad()
{

}

int main() { // compare to which value?
	double a = 0;
	double b = 1;
	
	double y = simpson(function, 100000, 0, 1);
	std::cout << "y1 = " << std::setprecision(15) << y << " y2 = "  << factorial(10) << std::endl;
	printVector(nthlegendre(3));
}
