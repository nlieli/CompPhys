#include <iostream>
#include <complex>
#include "matplot/matplot.h"
#include "fftw3.h"
#include "SF2.h"
#include "VectorMath.h"

#define EXERCISE 13

namespace ct
{
	const double PI =  3.1415926535897932;
	const std::string fileName = "txt_sound_file.txt";
}

static std::vector<std::complex<double>> discreteFourierTransform(std::vector<double> function, Timer* timer = nullptr)
{
	if (timer)
		timer->start();

	int N = (int)function.size();
	std::vector<std::complex<double>> result(N);
	std::complex<double> i(0.0, 1.0);

	for (size_t k = 0; k < N; ++k)
	{
		for (size_t j = 0; j < N; ++j)
			result[k] += function[j] * std::exp(-i * (std::complex<double>)(2 * ct::PI * j * k / N));
	}

	if (timer)
		timer->end();

	return result;
}

static std::vector<double> powerSpectrum(std::vector<std::complex<double>> function)
{
	size_t N = function.size();
	std::vector<double> result(N);
	double N2 = (double)N * N;
	for (size_t i = 0; i < N; ++i)
		result[i] = std::norm(function[i]) / N;

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
		std::cout << "\033[1;35m---------- Exercise 1a) ----------\033[0m" << std::endl;
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
		Timer timerDFT;
		Timer timerFFT;
		using matrix = std::vector<std::vector<double>>;
		int m = 5e2;
		matrix data = readMatrix(ct::fileName);
		std::vector<std::complex<double>> dft;
		std::vector<std::complex<double>> fft;
		for (int i = m; i > 0; --i)
		{
			data[0].resize(i);
			dft.resize(i);
			dft = discreteFourierTransform(data[0], &timerDFT);
			fft = fftw3::fft(data[0], &timerFFT);
		}

		std::vector<float> dftTimer = timerDFT.laps;
		std::vector<float> fftTimer = timerFFT.laps;
		std::reverse(dftTimer.begin(), dftTimer.end());
		std::reverse(fftTimer.begin(), fftTimer.end());

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
			grid(on);
			show();
		}
#endif
	}
#endif

// ------ 1c & d) ------
#if EXERCISE == 13
	{
		using matrix = std::vector<std::vector<double>>;
		matrix data = readMatrix(ct::fileName);
		int samplingRate = 44100;
		size_t N = data[0].size();

		std::vector<double> frequency = linspace(1, N, 0);
		frequency = frequency * samplingRate / N;

		std::vector<std::complex<double>> fft = fftw3::fft(data[0]);
		std::vector<double> PS = powerSpectrum(fft);

		std::vector<double>::iterator max = std::max_element(PS.begin(), PS.end());
		int dist = std::distance(PS.begin(), max);

		std::cout << "\033[1;35m---------- Exercise 1c & d) ----------\033[0m" << std::endl;
		std::cout << "The note played by the sound file is \033[1;31mD3\033[0m and its fundamental frequency is 146.83 Hz.\n"
			<< "The following result is given for the fundamental frequency by the FFT "
			<< "f = \033[1;31m" << frequency[dist] << " Hz\033[0m" << std::endl;


#ifdef NDEBUG
		{
			using namespace matplot;
			figure();
			plot(frequency, PS);
			xlabel("f / Hz");
			ylabel("Magnitude / 1");
			grid(on);
			show();
		}

#endif

	}
#endif

}
