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

double SINGULARITY_ERROR_TERM = 1e-5; // smaller term worsens error for small N but can increase accuracy for large N
namespace ct
{
	const double pi = 3.1415926535897932;
}

// ---------- 1a) specific functions ----------
static double F1a_integrand(double x)
{
	return 1 / std::sqrt(1 - std::pow( x, 4));
}

// ---------- 1b) specific functions ----------
static double T1b_integrand(double(*potential)(double x), double a, double x)
{
	return 1 / std::sqrt(potential(a) - potential(x));
}

static double V1(double x) { return std::cosh(x); }

static double V2(double x) { return std::exp(std::abs(x)); }

static double V3(double x) { return -std::cos(x); }

// ---------- 1c) specific functions ----------
static double T1c_integrand(double(*potential)(double x), double a, double k, double x)
{
	return 1 / std::sqrt(potential(a) - potential(k * x));
}

static double Vt(double x) { return std::tanh(x); }

// ---------- general purpose functions ----------
static double trapezoid(std::function<double(double)> integrand, int N, double a, double b)
{
	if (isinf(integrand(a))) a += SINGULARITY_ERROR_TERM; 
	if (isinf(integrand(b))) b -= SINGULARITY_ERROR_TERM;
	long double dx = (b - a) / N;
	long double F = 0.5 * (integrand(a) + integrand(b));

	for (int i = 1; i < N; ++i)
	{
		F += integrand(a + i * dx);
	}

	return F * dx;
}

static double simpson(std::function<double(double)> integrand, int n, double a, double b)
{
	if (isinf(integrand(a))) a += SINGULARITY_ERROR_TERM;
	if (isinf(integrand(b))) b -= SINGULARITY_ERROR_TERM;
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

	f = dx / 3 * (integrand(a0) + integrand(a + dx) + 4 * S1 + 2 * S2);
	
	return f;
}

static double legendreGauss(std::function<double(double)> integrand, int N, double a, double b)
{
	double F = 0;
	double chi;
	int RIT = 10; // number of Newton Raphson iterations
	std::vector<double> w(N);	
	std::vector<double> x(N);
	
	x = legendreFindRoots(N, RIT, -1, 1);
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
{ 
	// ---------- 1a) ----------
	{ /*
		double a = 0; // lower bound for integration
		double b = 1; // upper bound for integration 
		int N = 100; // number of integration range subdivisions

		std::vector<std::vector<double>> results(3, std::vector<double>(N));

		for (int i = 1; i < N + 1; ++i)
		{
			results[0][i - 1] = trapezoid(F1a_integrand, i, a, b);
			results[2][i - 1] = legendreGauss(F1a_integrand, i, a, b);
		}

		int j = 0;
		for (int i = 1; i < 2 * N; i += 2)
		{
			results[1][j] = simpson(F1a_integrand, i, a, b);
			j += 1;
		}

		plotResults(results);
		*/
	}
	// ---------- 1b) ----------
	{ /*
		std::vector<double> a = linspace(0, ct::pi, 100);
		size_t length = a.size();
		std::vector<double> Ta1(length);
		std::vector<double> Ta2(length);
		std::vector<double> Ta3(length);
	
		for (int i = 0; i < length; ++i)
		{
			double ap = a[i];
			Ta1[i] = legendreGauss([ap](double x) {return T1b_integrand(V1, ap, x); }, 10, 0, ap);
			Ta2[i] = legendreGauss([ap](double x) {return T1b_integrand(V2, ap, x); }, 10, 0, ap);
			Ta3[i] = legendreGauss([ap](double x) {return T1b_integrand(V3, ap, x); }, 10, 0, ap);
		}

		matplot::figure();
		matplot::plot(a, Ta1);
		matplot::hold(matplot::on);
		matplot::plot(a, Ta2);
		matplot::plot(a, Ta3);
		matplot::legend({ "cosh(x)", "exp(|x|)", "-cos(x)" });
		matplot::show();
		*/	
	}
	// ---------- 1c) ----------
	{
		double a = 1;
		double kp;
		std::vector<double> k = linspace(-20, 0, 100);
		size_t length = k.size();
		std::vector<double> Ta(length);
		std::vector<double> Ta2(length);

		for (int i = 0; i < length; ++i)
		{
			kp = k[i];
			Ta[i] = legendreGauss([a, kp](double x) {return T1c_integrand(Vt, a, kp, x); }, 10, 0, a);
			Ta2[i] = legendreGauss([a, kp](double x) {return T1c_integrand(Vt, a / 2, kp, x); }, 10, 0, a / 2);
		}
		
		printVector(k);
		printVector(Ta);
		matplot::figure();
		matplot::plot(k, Ta);
		matplot::hold(matplot::on);
		matplot::plot(k, Ta2);
		matplot::show();
	}


	/*
	double mu = 3e-6;
	std::vector<double> f = {-mu, 2 * mu, -mu, (3 - 2 * mu), (mu - 3), 1};
	std::vector<double> x = matplot::linspace(-1, 2);
	std::vector<double> y(x.size());
	std::vector<double> roots;

	for (int i = 0; i < x.size(); ++i) 
	{
		y[i] = polyEval(f, x[i]);
	}

	roots = polyFindRoots(f, 10, -1, 1);
	printVector(roots);
	matplot::figure();
	matplot::plot(x, y);
	matplot::axis(matplot::equal);
	matplot::grid(matplot::on);
	matplot::show();
	*/
	

	//printVector(y);

	//printVector(results[0]);
	//printVector(results[1]);
	//printVector(results[2]);

	//T_integrand(cosh, 2, 3);
	//T_integrand(exp, abs(a), abs(3));
	//T_integrand([] (double x) {return -cos(x); }, a, 3);


	
}
