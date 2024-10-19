#include "recursive.h"
#include <vector>
#include <immintrin.h>
#include "SF.h"

double dx1 = 1e-5;

std::vector<double> legendreBracketing(int N, double a, double b)
{
	std::vector<double> guess;
	double i = a;
	for (double j = a; j < b; j += dx1)
	{
		if (legendreBrackSearch(N, i, j))
		{
			guess.push_back(i + dx1 / 2);
		}
		i = j;
	}
	
	return guess;
}


bool legendreBrackSearch(int N, double a, double b)
{
	double y1 = legendreEval(N, a);
	double y2 = legendreEval(N, b);
	double epsilon = 1e-10;

	if ((y1 * y2) < 0)
		return true;

	if (std::abs(y1) < epsilon || std::abs(y2) < epsilon)
		return true;

	return false;
}

std::vector<double> legendreFindRoots(int N, int RIT, double a, double b)
{
	std::vector<double> guess = legendreBracketing(N, a, b);
	size_t numberOfRoots = guess.size();
	std::vector<double> result(numberOfRoots);
	double epsilon = 1e-100;
	double proximityCheck = 0;

	for (size_t i = 0; i < numberOfRoots; ++i)
		result[i] = legendreNewtonRaphson(N, RIT, guess[i]);

	for (size_t i = 0; i < numberOfRoots - 1; ++i)
	{
		proximityCheck = result[i] - result[i + 1];
		if (std::abs(proximityCheck) < epsilon)
		{
			result.erase(result.begin() + i);
			numberOfRoots -= 1;
			i -= 1;
		}
	}

	return result;
}

double legendreNewtonRaphson(int N, int RIT, double x0)
{
	for (int i = 0; i < RIT; ++i)
	{
		if (x0 == 0)
			return x0;
		x0 = x0 - (legendreEval(N, x0) / legendrePrimeEval(N, x0));
	}

	return x0;
}

std::vector<double> legendreWeights(std::vector<double> xi, int N)
{
	std::vector<double> weights(N);

	for (int i = 0; i < N; ++i)
	{
		double PnPrime = legendrePrimeEval(N, xi[i]);
		weights[i] = 2 / ((1 - xi[i] * xi[i]) * PnPrime * PnPrime);
	}

	return weights;
}




