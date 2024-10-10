#include <vector>
#include <iostream>
#include <cmath>

unsigned int factorial(unsigned int x) // replace with gamma function ?? 
{
	int y = 1;

	for (int i = x; i > 1; --i)
		y *= i;

	return y;
}

unsigned int nCr(unsigned int n, unsigned int r) // generalize with templates to take size_t
{
	return factorial(n) / (factorial(r) * factorial(n - r));
}

std::vector<double> nthlegendre(int n)
{
	std::vector<double> Pn(n + 1);
	for (int k = 0; k <= n / 2; ++k)
		Pn[n - 2 * k] = std::pow(-1, k) * nCr(n, k) * nCr(2 * n - 2 * k, n) * std::pow(2, -n);
	return Pn;
}

