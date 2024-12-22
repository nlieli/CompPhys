#include <iostream>
#include <complex>
#include "matplot/matplot.h"
#include "fftw3.h"
#include "SF2.h"
#include "VectorMath.h"

#define EXERCISE 12

namespace ct
{
	const double PI =  3.1415926535897932;
	const std::string fileName = "txt_sound_file.txt";
}

static std::vector<std::complex<double>> discreteFourierTransform(std::vector<double> function, Timer* timer = nullptr)
{
	if (timer)
		timer->start();

	std::vector<std::complex<double>> result(function.size());
	std::complex<double> i(0.0, 1.0);
	int N = (int)function.size();

	for (size_t k = 0; k < N; ++k)
	{
		for (size_t j = 0; j < N; ++j)
			result[k] += function[j] * std::exp(-i * (std::complex<double>)(2 * ct::PI * j * k / N));
	}

	if (timer)
		timer->end();

	return result;
}


int main()
{

// ------ 1a) ------
#if EXERCISE == 11
	{
		using matrix = std::vector<std::vector<double>>;
		matrix data = readMatrix(ct::fileName);

		std::vector<std::complex<double>> dft = discreteFourierTransform(data[0]);
		std::vector<std::complex<double>> fft = fftw3::fft(data[0]);
		std::vector<std::complex<double>> diff = dft - fft;

		std::complex<double> sum = 0;
		double epsilon = 1e-8;
		for (size_t i = 0; i < diff.size(); ++i)
		{
			double u = (std::abs(diff[i].real()) < epsilon) ? 0.0 : diff[i].real();
			double v = (std::abs(diff[i].imag()) < epsilon) ? 0.0 : diff[i].imag();

			diff[i] = std::complex<double>(u, v);
			sum += diff[i];
		}
		std::cout << "\033[1;35m---------- Exercice 1a) ----------\033[0m" << std::endl;
		std::cout << "Sum of differences between DFT and FFT = "
			<< "\033[3;36m"
			<< sum.real()
			<< " + i"
			<< sum.imag()
			<< "\033[0m"
			<< " using epsilon = "
			<< "\033[3;36m"
			<< epsilon
			<< "\033[0m"
			<< std::endl;


#ifdef NDEBUG
		{
			using namespace matplot;
			figure();
			plot(data[0]);
			hold(on);
			plot(data[1]);

			show();
		}
#endif
	}
#endif

// ------ 1b) ------
#if EXERCISE == 12
	{
		Timer timer;
		using matrix = std::vector<std::vector<double>>;
		int m = 5e2;
		matrix data = readMatrix(ct::fileName);
		std::vector<std::complex<double>> dft;
		std::vector<std::complex<double>> fft;
		std::vector<double> dataVec = data[0];
		for (int i = m; i > 0; --i)
		{
			data[0].resize(i);
			dft.resize(i);
			dft = discreteFourierTransform(data[0], &timer);
		}

		std::vector<float> dftTimer = timer.laps;
		std::reverse(dftTimer.begin(), dftTimer.end());

		for (int i = m; i > 0; --i)
		{
			dataVec.resize(i);
			fft.resize(i);
			fft = fftw3::fft(dataVec, &timer);
		}

		std::vector<float> fftTimer = timer.laps;
		std::reverse(fftTimer.begin(), fftTimer.end());
		fftTimer.resize(m);

		std::vector<float> dftTimeSum = cumsumVector(dftTimer);
		std::vector<float> fftTimeSum = cumsumVector(fftTimer);

#ifdef NDEBUG
		{
			using namespace matplot;
			figure();
			plot(dftTimeSum);
			hold(on);
			plot(fftTimeSum);
			xlabel("m");
			ylabel("ms");
			show();
		}
#endif
	}
#endif

// ------ 1c) ------
#if EXERCISE == 13
	{

#ifdef NDEBUG

#endif

	}
#endif

}