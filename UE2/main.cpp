#include <iostream>
#include <complex>
#include "matplot/matplot.h"
#include "fftw3.h"
#include "SF2.h"

#define EXERCICE 0

namespace ct
{
	const double PI =  3.1415926535897932;
}


static std::vector<std::complex<double>> discreteFourierTransform(std::vector<double> function)
{
	std::vector<std::complex<double>> result(function.size());
	std::complex<double> i(0.0, 1.0);
	int N = function.size();

	for (size_t k = 0; k < N; ++k)
	{
		for (size_t j = 0; j < N; ++j)
			result[k] += function[j] * std::exp(-i * (std::complex<double>)(2 * ct::PI * j * k / N));
	}

	return result;
}


int main()
{
#if EXERCICE == 0
	Timer timer;
	using matrix = std::vector<std::vector<double>>;
	matrix data = readMatrix("txt_sound_file.txt");
	std::vector<double> time = linspace(0, data[0].size(), 0);
	std::vector<double> test = { 1,2,3,4 };

	std::vector<std::complex<double>> dft = discreteFourierTransform(test);
	std::vector<std::complex<double>> fft = fftw3::fft(test);
	printVector(fft);

#ifdef NDEBUG
	{
		using namespace matplot;
		figure();
		plot(time, data[0]);
		hold(on);
		plot(time, data[1]);

		show();
	}
#endif
#endif

}
