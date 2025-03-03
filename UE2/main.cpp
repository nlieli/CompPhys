#include <iostream>
#include <complex>
#include "matplot/matplot.h"
#include "fftw3.h"
#include "SF2.h"
//#include "VectorMath.h"
#include <vector>
#include <algorithm>
#include "nstd.h"

#define nstd_print(var) nstd::print(var, #var) // calls print with variable name in output
#define EXERCISE 13 

namespace ct
{
	const double PI = 3.1415926535897932;
	const std::string fileName = "txt_sound_file.txt";
	const std::string fileName2 = "dna_coordinates.txt";
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

template <typename T, size_t S>
struct TriDiSoQ 
{
	std::array<std::array<T, S>, S> m_matrix = { 0 };
    std::array<T, S> m_rhsVec;
	std::array<T, S> u;

    TriDiSoQ(const std::array<T, S + 2>& x, const std::array<T, S + 2>& y) {
        std::array<T, S + 1> h;
        std::array<T, S> v;
        std::array<T, S + 1> b;
        //std::array<T, S> u;

        for (size_t i = 0; i < S + 1; ++i)
            h[i] = x[i + 1] - x[i];

        for (size_t i = 0; i < S; ++i)
            v[i] = 2 * (h[i] + h[i + 1]);

        for (size_t i = 0; i < S + 1; ++i)
            b[i] = (y[i + 1] - y[i]) / h[i];

        for (size_t i = 0; i < S; ++i)
            u[i] = 6 * (b[i + 1] - b[i]);

        nstd::diagonal_iterator it(m_matrix);
        nstd::diagonal_iterator it1(m_matrix, 1);
        nstd::diagonal_iterator it2(m_matrix, -1);

        for (size_t i = 0; it != it.end(); ++i) {
            *it = v[i]; 
            ++it;
        }

        for (size_t i = 0; it1 != it1.end(); ++i) {
            *it1 = h[i];  
            *it2 = *it1;  
            ++it1;
            ++it2;
        }

        m_rhsVec = u;
    }

	std::array<T, S + 2> solveLUDecomposition()
	{
		std::array<T, S> y = { 0 };
		std::array<T, S + 2> z = { 0 };
		nstd::diagonal_iterator main_it(m_matrix);
		nstd::diagonal_iterator super_it(m_matrix, 1);
		nstd::diagonal_iterator sub_it(m_matrix, -1);

		std::array<T, S> main_diag;
		std::array<T, S - 1> super_diag;
		std::array<T, S - 1> sub_diag;
		
		for (size_t i = 0; main_it != main_it.end(); ++main_it, ++i)
			main_diag[i] = *main_it;

		for (size_t i = 0; super_it != super_it.end(); ++super_it, ++i)
			super_diag[i] = *super_it;

		for (size_t i = 0; sub_it != sub_it.end(); ++sub_it, ++i)
			sub_diag[i] = *sub_it;

		std::array<T, S> upper_main_diag;
		std::array<T, S - 1> upper_super_diag;
		std::array<T, S - 1> lower_sub_diag;

		upper_main_diag[0] = main_diag[0];
		upper_super_diag = super_diag;
		
		for (size_t i = 1; i < S; ++i)
		{
			lower_sub_diag[i - 1] = upper_super_diag[i - 1] / upper_main_diag[i - 1];
			upper_main_diag[i] = main_diag[i] - upper_super_diag[i - 1] * lower_sub_diag[i - 1];
		}
		
		y[0] = m_rhsVec[0];
		for (size_t i = 1; i < S; ++i)
		{
			y[i] = m_rhsVec[i] - lower_sub_diag[i - 1] * y[i - 1];
		}

		z[S] = y[S - 1] / upper_main_diag[S - 1];
		for (int i = S - 2; i >= 0; --i)
		{
			z[i + 1] = (y[i] - upper_super_diag[i] * z[i + 2]) / upper_main_diag[i];
		}

		return z;
	}
};

template <typename T, size_t S>
static spline createSplines(const std::array<T, S>& x, const std::array<T, S>& y, const std::array<T, S>& z, size_t N)
{
	spline results;
	std::vector<double> ysegments;
	std::vector<double> xsegments;
	std::array<T, S - 1> h;
	for (size_t i = 0; i < S - 1; ++i)
		h[i] = x[i + 1] - x[i];

	for (size_t i = 0; i < S - 2; ++i)
	{
		auto xValues = nstd::linspace<double>(x[i], x[i + 1], N);
		ysegments = -z[i] / (6 * h[i]) * nstd::elemPowNdimArray((x[i + 1] - xValues), 3) +
		     z[i + 1] / (6 * h[i]) * nstd::elemPowNdimArray((xValues - x[i]), 3) +
			(y[i + 1] / h[i] - z[i + 1] / 6 * h[i]) * (xValues - x[i]) -
			(y[i] / h[i] - h[i] / 6 * z[i]) * (x[i + 1] - xValues);
		xsegments = xValues;
		results.yValues.insert(results.yValues.end(), ysegments.begin(), ysegments.end());
		results.xValues.insert(results.xValues.end(), xsegments.begin(), xsegments.end());
	}

	return results;
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
			xlabel("f / Hz");
			ylabel("Amplitude / 1");
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

		std::vector<float> dftTimeSum = nstd::cumsumVector(dftTimer);
		std::vector<float> fftTimeSum = nstd::cumsumVector(fftTimer);

#ifdef NDEBUG
		{
			using namespace matplot;
			figure();
			plot(dftTimeSum);
			hold(on);
			plot(fftTimeSum);
			xlabel("m");
			ylabel("ms");
			legend("DFT", "FFT");
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

		std::vector<double> frequency = nstd::linspace(1.0, (double)N, 0);
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
			hold(on);
			
			xlabel("f / Hz");
			ylabel("Magnitude / 1");
			grid(on);
			show();
		}

#endif

	}
#endif

	// ------ 1e) ------
#if EXERCISE == 14
	{
		std::cout << "\033[1;35m---------- Exercise 1e) ----------\033[0m" << std::endl;
		std::cout << "The Fourier representation of a perfect tone would be a Dirac delta function delta(x - f_k) with f_k the fundamental frequency of the tone."
			<< std::endl;

		using matrix = std::vector<std::vector<double>>;
		constexpr size_t numberOfPeaks = 5;
		matrix data = readMatrix(ct::fileName);
		std::vector<std::complex<double>> approximatedFFT(data[0].size(), 0); // creates a vector with length data[0].size() filled with all zeros

		std::vector<std::complex<double>> fft = fftw3::fft(data[0]);
		std::vector<double> PSe = powerSpectrum(fft);
		PSe.resize(PSe.size() / 2); // removes values above Nyquist Frequency containing no additional information

		// find the 5 tallest peaks and copies the info to the approximation (contains phase information as well)
		int approximatePeakWidth = 30; // catches all data points belonging to one peak
		for (size_t i = 0; i < numberOfPeaks; ++i)
		{
			std::vector<double>::iterator max = std::max_element(PSe.begin(), PSe.end());
			int index = std::distance(PSe.begin(), max);
			approximatedFFT[index] = fft[index];
			std::fill(PSe.begin() + index - approximatePeakWidth, PSe.begin() + index + approximatePeakWidth, 0); // removes largest peak from data set so second largest can be found
		}

		std::vector<double> approxmiatedSoundWave = fftw3::ifft(approximatedFFT);

		/* use matlab script to actually play the sound audio
		- i cant be bother to install another library and figure
		out how it works just for that

		the following Matlab code will give you a correct output

		data = readmatrix("approximatedFFT.txt")
		Fs = 44100;
		sound(data, Fs);

		just copy and paste, then you will hear that the tone is a
		perfect D3.

		alternatively, look at the wave form of the file in the
		graph.
		*/


		// writes output file with sound wave values
		std::ofstream fs;
		std::string fileName = "approximatedFFT.txt";
		fs.open(fileName);

		for (size_t i = 0; i < approxmiatedSoundWave.size(); ++i)
			fs << approxmiatedSoundWave[i] << "\n";

		fs << std::endl;
		fs.close();

#ifdef NDEBUG
		{
			using namespace matplot;
			figure();
			plot(approxmiatedSoundWave);
			xlabel("t");
			ylabel("Amplitude");
			show();
		}
#endif
	}
#endif

	// ------ 2) ------
#if EXERCISE == 2
	{
		const size_t xy_i = 11;
		const size_t hv_i = 9;
		using matrix = std::array<std::array<double, hv_i>, hv_i>;
		std::array<double, xy_i> x = { -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5 };
		std::array<double, xy_i> y;
		for (size_t j = 0; j < xy_i; ++j)
			y[j] = 1 / (1 + (double)x[j] * x[j]);

		TriDiSoQ<double, hv_i> system(x, y);
		std::array<double, xy_i> z = system.solveLUDecomposition();
		spline splines = createSplines(x, y, z, 100);

		size_t relaxations_values = 100;
		auto omega = nstd::linspace(0.1, 1.9, relaxations_values);

		matrix A = system.m_matrix;
		for (nstd::diagonal_iterator it(A); it != it.end(); ++it)
			*it = 0;

		std::vector<double> iterations(relaxations_values, 0.0);
		std::array<double, xy_i> z_approximation = { 0 };
		std::array<double, hv_i> z_approximation_small = { 0 };
		std::array<double, xy_i> diff;
		double tolerance = 1e-12;

		for (size_t i = 0; i < relaxations_values; ++i)
		{
			do
			{
				for (size_t j = 0; j < hv_i; ++j)
				{
					z_approximation[j + 1] = (1 - omega[i]) * z_approximation[j + 1] + omega[i] / system.m_matrix[j][j] *
						(system.u[j] - A[j] * z_approximation_small);
					std::copy(z_approximation.begin() + 1, z_approximation.end() - 1, z_approximation_small.begin());
				}
				++iterations[i];
				diff = nstd::absNdimArray(z_approximation - z);
			} while (*std::max_element(diff.begin(), diff.end()) > tolerance);
			std::fill(z_approximation.begin(), z_approximation.end(), 0);
			std::fill(z_approximation_small.begin(), z_approximation_small.end(), 0);
		}

#ifdef NDEBUG
		{
			using namespace matplot;
			figure();
			plot(splines.xValues, splines.yValues);
			figure();
			plot(iterations);
			show();
		}
#endif
	}
#endif

	// ------ 3) ------
#if EXERCISE == 3
	{
		using matrix = std::vector<std::vector<double>>;
		matrix data = readMatrix(ct::fileName2, comma);
		std::vector<double> x = data[0];
		std::vector<double> y = data[1];
		std::vector<double> z = data[2];
		std::vector<double> m = data[3];
		size_t N = data[0].size();
		matrix H(N, std::vector<double>(N, 0));
		matrix M(N, std::vector<double>(N, 0));
		matrix mInvRoots(N, std::vector<double>(N, 0));

		const int R_cut= 5;
		const int k = 1;

		double dx, dy, dz;
		double R;
		int neighbouringAtoms = 0;

		for (size_t i = 0; i < N; ++i)
		{
			neighbouringAtoms = 0;
			for (size_t j = 0; j < N; ++j)
			{
				if (i != j)
				{
					dx = x[i] - x[j];
					dy = y[i] - y[j];
					dz = z[i] - z[j];

					R = std::sqrt(dx * dx + dy * dy + dz * dz);
					if (R < R_cut)
					{
						++neighbouringAtoms;
						H[i][j] = -k;
					}
				}
			}

			H[i][i] = neighbouringAtoms * k;
			M[i][i] = m[i];
			mInvRoots[i][i] = 1 / std::sqrt(m[i]);
		}
	
		size_t numberOfEigenVectors = 10;
		matrix K = mInvRoots * H * mInvRoots;
		matrix K_deflation(N, std::vector<double>(N));
		matrix eigenVectors(numberOfEigenVectors, std::vector<double>(N));

		size_t iterations = 1000;
		std::vector<double> lambda(numberOfEigenVectors);
		std::vector<double> powerMethodVec(N, 1.0);

		for (size_t i = 0; i < numberOfEigenVectors; ++i)
		{
			for (size_t j = 0; j < iterations; ++j)
			{
				powerMethodVec = K * powerMethodVec;
				powerMethodVec = 1 / std::sqrt(powerMethodVec * powerMethodVec) * powerMethodVec;
			}
			lambda[i] = powerMethodVec * (K * powerMethodVec) / (powerMethodVec * powerMethodVec);
			eigenVectors[i] = powerMethodVec;
			K = K - lambda[i] * nstd::vectorProdToMatrix(powerMethodVec, powerMethodVec);
		}
		matrix firstStrand(eigenVectors.size(), std::vector<double>(eigenVectors[0].size() / 2));
		matrix secondStrand(eigenVectors.size(), std::vector<double>(eigenVectors[0].size() / 2 + 1));
		std::vector<double> zFirstStrand(z.size() / 2);
		std::vector<double> zSecondStrand(z.size() / 2 + 1);

		for (size_t i = 0; i < eigenVectors.size(); ++i) {
			std::copy(eigenVectors[i].begin(), eigenVectors[i].begin() + eigenVectors[i].size() / 2, firstStrand[i].begin());
			std::copy(eigenVectors[i].begin() + eigenVectors[i].size() / 2, eigenVectors[i].end(), secondStrand[i].begin());
		}

		std::copy(z.begin(), z.begin() + z.size() / 2, zFirstStrand.begin());
		std::copy(z.begin() + z.size() / 2, z.end(), zSecondStrand.begin());

#ifdef NDEBUG
		{
			using namespace matplot;
			std::cout << "\033[1;35m---------- Exercise 3) ----------\033[0m" << std::endl;
			std::cout << "First Strand is \033[1;38;2;0;150;255mblue\033[0m" << std::endl;
			std::cout << "Second Strand is \033[1;38;2;255;115;0morange\033[0m" << std::endl;

			figure();
			for (size_t i = 0; i < 4; ++i)
			{
				subplot(2, 2, i);
				plot(zFirstStrand, firstStrand[i]);
				hold(on);
				plot(zSecondStrand, secondStrand[i]);
				std::string stitle = "lambda = ";
				std::string lm = std::to_string(lambda[i]);
				stitle = stitle + lm;
				title(stitle);
			}

			figure();
			for (size_t i = 4; i < 8; ++i)
			{
				subplot(2, 2, i);
				plot(zFirstStrand, firstStrand[i]);
				hold(on);
				plot(zSecondStrand, secondStrand[i]);
				std::string stitle = "lambda = ";
				std::string lm = std::to_string(lambda[i]);
				stitle = stitle + lm;
				title(stitle);
			}

			figure();
			for (size_t i = 8; i < 10; ++i)
			{
				subplot(1, 2, i);
				plot(zFirstStrand, firstStrand[i]);
				hold(on);
				plot(zSecondStrand, secondStrand[i]);
				std::string stitle = "lambda = ";
				std::string lm = std::to_string(lambda[i]);
				stitle = stitle + lm;
				title(stitle);
			}
			show();
		}
#endif
	}
#endif

}
