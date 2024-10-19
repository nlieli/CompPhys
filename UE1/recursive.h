#pragma once
#include <vector>


std::vector<double> legendreBracketing(int N, double a, double b);

bool legendreBrackSearch(int N, double a, double b);

std::vector<double> legendreFindRoots(int N, int RIT, double a, double b);

double legendreNewtonRaphson(int N, int RIT, double x0);

std::vector<double> legendreWeights(std::vector<double> xi, int N);



template <typename T>
T legendreEval(int N, T x)
{
	T Pn_ = 1;
	T Pn = x;
	T Pnx = 0;

	if (N == 0)
		return Pn_;

	if (N == 1)
		return Pn;
	
	for (int i = 1; i < N; ++i)
	{
		Pnx = ((2 * i + 1) * x * Pn - i * Pn_) / (i + 1);
		Pn_ = Pn;
		Pn = Pnx;
	}

	return Pn;
}

template <typename T>
T legendrePrimeEval(int N, T x)
{
	T Pn = legendreEval(N, x);
	T Pn_ = legendreEval(N - 1, x);

	return N / (x * x - 1) * (x * Pn - Pn_);
}


