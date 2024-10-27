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
#include <random> // pseudo random device

// Timer in debug mode if needed
#ifdef _DEBUG
#define TIMER Timer timer
#else 
#define TIMER
#endif

/* --- !IMPORTANT! ---
To select the specific excercie you want to look at, you need to change the definition of the
EXERCIC macro variable. You may change the number to one of the following:
#define EXERCICE 0 --- corresponds to viewing all exercices
#define EXERCICE 11 --- corresponds to viewing 1a)

all possible values: 0 = all, 11 = 1a), 12 = 1b), 13 = 1c), 2 = 2a) and 2b), 3 = 3a) and 3b)
When using Visual Studio, all section that are not selected appear as greyed out. If you define it
as an invalid value, the preprocessor will remove all the code and none of the exercices will
give any output.
*/

#define EXERCICE 3

double SINGULARITY_ERROR_TERM = 1e-5; // smaller term worsens error for small N but can increase accuracy for large N
namespace ct
{
	const double PI = 3.1415926535897932;
	const double EARTH_MASS_KG = 5.97219e24;
	const double SUN_MASS_KG = 1.9891e30;
	const double DISTANCE_EARTH_SUN_METERS = 149e6;
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
	return 1 / std::sqrt(potential(k * a) - potential(k * x));
}

static double Vt(double x) { return (std::tanh(x)); }

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

static double simpson(std::function<double(double)> integrand, int N, double a, double b)
{
	if (isinf(integrand(a))) a += SINGULARITY_ERROR_TERM;
	if (isinf(integrand(b))) b -= SINGULARITY_ERROR_TERM;
	double dx = (b - a) / N;
	double f = 0;
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
#if EXERCICE == 11 || EXERCICE == 0
	{ 
		double a = 0; // lower bound for integration
		double b = 1; // upper bound for integration 
		int N = 100; // number of integration range subdivisions
		double trueValue = 1.311028777146120;

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

		//plotResults(results);
		
		results[0] = abs(results[0] - trueValue); // trap
		results[1] = abs(results[1] - trueValue); // simp
		results[2] = abs(results[2] - trueValue); // gauss

		std::vector<double> GTiterations = linspace(0, results[0].size(), 0); // 0 means step width of int 1
		std::vector<double> simpIterations = linspace(0, results[1].size(), 0); // simpson is different lenght because only even Ns allowed

#ifdef NDEBUG
		{
			using namespace matplot;
			semilogy(GTiterations, results[0]);
			hold(on);
			semilogy(simpIterations, results[1]);
			semilogy(GTiterations, results[2]);
			auto lg = matplot::legend({ "Trapezoid", "Simpson", "Gaussian Quadrature" });
			lg->font_name("Arial");
			title("1a)");
#if EXERCICE == 11
			show();
#endif
		}
#endif
	}
#endif
	// ---------- 1b) ----------
#if EXERCICE == 12 || EXERCICE == 0
	{ 
		std::vector<double> a = linspace(0, ct::PI, 100);
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

#ifdef NDEBUG
		{
			using namespace matplot;
			figure();
			plot(a, Ta1);
			hold(matplot::on);
			plot(a, Ta2);
			plot(a, Ta3);
			auto lg = matplot::legend({ "cosh(x)", "exp(|x|)", "-cos(x)" });
			lg->font_name("Arial");
			title("1b");
#if EXERCICE == 12
			show();
#endif
		}
#endif	
	}
#endif
	// ---------- 1c) ----------
#if EXERCICE == 13 || EXERCICE == 0
	{ 
		double a = 1;
		double kp;
		std::vector<double> k = linspace(0, 10, 100);
		size_t length = k.size();
		std::vector<double> Ta(length);
		std::vector<double> Ta2(length);
		std::vector<double> TaR(length);

		for (int i = 0; i < length; ++i)
		{
			kp = k[i];
			Ta[i] = legendreGauss([a, kp](double x) {return T1c_integrand(Vt, a, kp, x); }, 10, 0, a);
			Ta2[i] = legendreGauss([a, kp](double x) {return T1c_integrand(Vt, a / 2, kp, x); }, 10, 0, a / 2);
			TaR[i] = Ta[i] / Ta2[i];
		}
		
		//printVector(k);
		//printVector(Ta);

#ifdef NDEBUG
		{
			using namespace matplot;
			figure();
			plot(k, Ta);
			hold(matplot::on);
			plot(k, Ta2);
			title("1c) 1");

			figure();
			plot(k, TaR);
			title("1c) 2");
#if EXERCICE == 13
			show();
#endif
		}
#endif
	}
#endif
	// ---------- 2) ----------
#if EXERCICE == 2 || EXERCICE == 0
	{
		const int iterations = 16;
		const double trueValue = 0.009969265283747;
		const double mu = ct::EARTH_MASS_KG / (ct::EARTH_MASS_KG + ct::SUN_MASS_KG);
		std::vector<double> quintic = { -mu, 2 * mu, -mu, (3 - 2 * mu), (mu - 3), 1 };

		std::vector<std::array<double, 2>> interval = polyBracketingInterval(quintic, -10, 10); // vector in case there are more than one root
		//std::vector<std::array<double, 2>> interval = { {-10, 10} }; // adjust values to edit initial guess
		size_t numberOfRoots = interval.size();
		std::vector<std::array<double, iterations>> newtonResult(numberOfRoots);
		std::vector<std::array<double, iterations>> bisectionResult(numberOfRoots);
		std::array<double, iterations> bisectionRoot; // approximation of root for each iteration
		std::array<double, iterations> newtonRoots;

		for (size_t i = 0; i < numberOfRoots; ++i)
		{
			for (int j = 0; j < iterations; ++j)
			{
				newtonRoots[j] = polyNewtonRaphson(quintic, j, (interval[i][0] + interval[i][1]) / 2);
				bisectionRoot[j] = polyBisection(quintic, j, interval[i]);
			}

			newtonResult[i] = newtonRoots; // for each root, all approximations are stored in an array
			bisectionResult[i] = bisectionRoot;
		}
		double lagrangePoint = newtonResult[0][iterations - 1] * ct::DISTANCE_EARTH_SUN_METERS;
		std::cout << "Lagrange Point is at: " <<  lagrangePoint << " m from earth" << std::endl;
		
		std::cout << "Newton Root Convergence: " << std::endl;
		print(newtonResult[0] - trueValue);

		std::cout << "Bisection Root Convergence: " << std::endl;
		print(bisectionResult[0] - trueValue);

#ifdef NDEBUG
		{
			using namespace matplot;
			std::vector<double> n = linspace(0, iterations, 0);
			std::array<double, iterations> newton = arrayLog10(abs(newtonResult[0] - trueValue));
			std::array<double, iterations> bisection = arrayLog10(abs(bisectionResult[0] - trueValue));

			figure();
			plot(n, newton);
			hold(on);
			plot(n, bisection);
			auto lg = matplot::legend({ "Newton", "Bisection" });
			lg->font_name("Arial");
			xlabel("iteration");
			ylabel("Error / ln(|x - x_t|)");
			title("2a)");

			figure();
			plot(newton, arrayShiftLeft(newton)); // needs to be made more pretty
			title("2b)");
#if EXERCICE == 2
			show();
#endif
		}
#endif
	}
#endif
	// ---------- 3) ----------
#if EXERCICE == 3 || EXERCICE == 0
	{
		const int iterations = 1000;
		std::random_device rd;
		std::uniform_real_distribution<double> distribution(0, 1);
		std::vector<std::vector<double>> rootMatrix;
		rootMatrix.reserve(iterations);
		int numberOfRoots = 0;

		std::cout << "Loading... " << std::endl;
		for (int i = 0; i < iterations; ++i)
		{
			std::vector<double> polynomial(6);
			for (double& value : polynomial)
				value = distribution(rd);

			std::vector<double> roots = polyFindRoots(polynomial, 10, -10, 10);
			rootMatrix.emplace_back(roots);

			numberOfRoots += roots.size();
		}

		double meanRoots = (double)numberOfRoots / iterations;
		std::cout << "Average number of Roots in [-10, 10] = " << meanRoots << std::endl; // 3b)

		std::vector<double> allRoots;
		for (const std::vector<double>& rootVector : rootMatrix)
			allRoots.insert(allRoots.end(), rootVector.begin(), rootVector.end()); // values for 3a)


#ifdef NDEBUG
		{
		using namespace matplot;
		figure();
		hist(allRoots, 100);
		xlim({ -10, 10 });
		title("3a)");
		show();
		}
#endif
	}
#endif

}

