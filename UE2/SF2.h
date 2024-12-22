#pragma once
#include <fstream>

enum decimalStyle
{
	comma = ',',
	dot = '.'
};

namespace fftw3
{
	std::vector<std::complex<double>> fft(std::vector<double> values);
}
std::vector<std::vector<double>> readMatrix(const std::string& fileName, const decimalStyle& style = dot);

template <typename T> 
void printVector(std::vector<T>& vector)
{
	std::cout << '[';
	for (size_t i = 0; i < vector.size(); ++i)
		std::cout << vector[i] << ' ';
	std::cout << ']' << std::endl;
}

struct Timer
{
	std::chrono::time_point<std::chrono::steady_clock> start, end;
	std::chrono::duration<float> duration;

	Timer();

	~Timer();
};

template <typename T, typename B, typename C>
std::vector<double> linspace(T a, B b, C n = 0) // a = lower bound, b = higher bound, n = number of values
{
	if (n == 0) n = b - a;

	double increment = (double)(b - a) / n;
	std::vector<double> vector(n);

	for (int i = 0; i < n; ++i)
		vector[i] = a + i * increment;

	return vector;
}

