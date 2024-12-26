#include <iostream>
#include <fstream> // creates input file stream
#include <sstream>
#include <string>
#include <regex>
#include <chrono>
#include <vector>
#include <complex>
#include "fftw3.h"
#include "SF2.h"

Timer::Timer()
{
	global_start = std::chrono::high_resolution_clock::now();
	global_duration = std::chrono::duration<float> (0.0f);
	global_end = global_start;

	local_start = std::chrono::high_resolution_clock::now();
	local_end = local_start;
	local_duration = global_duration;
}

Timer::~Timer()
{
	global_end = std::chrono::high_resolution_clock::now();
	global_duration = global_end - global_start;

	float ms = global_duration.count() * 1000.0f;
	std::cout << "Total program timer took " << ms << "ms " << std::endl;
}

void Timer::start()
{
	local_start = std::chrono::high_resolution_clock::now();
}

void Timer::end()
{
	local_end = std::chrono::high_resolution_clock::now();
	local_duration = local_end - local_start;
	float ms = local_duration.count() * 1000.0f;
	laps.push_back(ms);
}

std::vector<std::vector<double>> readMatrix(const std::string& fileName, const decimalStyle& style)
{
	std::ifstream fs;
	std::stringstream ss;
	std::string row;
	std::vector<std::vector<double>> Matrix;
	std::regex reg;
	size_t columns = 0;

	if (style == dot)
		reg = "(-?\\d*.?\\d+e?[-+]\\d*)"; // decimal point system, captures all numbers using this style 
	else
		reg = "(-?\\d*,?\\d+e?[-+]\\d*)"; // decimal comma system, captures all numbers using this style

	fs.open(fileName, std::ios::in);
	if (!fs.is_open()) { std::cerr << "[ERROR] File stream did not open successfully!" << std::endl;  return Matrix; }

	ss << fs.rdbuf();
	std::getline(ss, row); // offsets ss pointer by one line
	std::smatch matches;
	std::regex_search(row, matches, reg);
	columns = matches.size();
	Matrix.resize(columns);
	double value;

	ss.clear();
	ss.seekg(0, std::ios::beg); // resets ss pointer to start of ss object to parse all the data

	for (int i = 0; ss >> value; ++i)
		Matrix[i % columns].push_back(value);
	
	fs.close();
	return Matrix;
}

std::vector<std::complex<double>> fftw3::fft(std::vector<double> values, Timer* timer)
{
	if (timer)
		timer->start();

	int N = (int)values.size();
	std::vector<std::complex<double>> results(N);
	fftw_complex* in = fftw_alloc_complex(N);
	fftw_complex* out = fftw_alloc_complex(N);

	for (int i = 0; i < N; ++i)
	{
		in[i][0] = values[i];
		in[i][1] = 0.0;
	}

	fftw_plan plan = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(plan);

	for (int i = 0; i < N; ++i)
		results[i] = std::complex<double>(out[i][0], out[i][1]);

	fftw_destroy_plan(plan);
	fftw_free(in);
	fftw_free(out);
	fftw_cleanup();

	if (timer)
		timer->end();

	return results;
}


std::vector<double> fftw3::ifft(std::vector<std::complex<double>> values)
{
	int N = (int)values.size();
	std::vector<double> results(N);
	fftw_complex* in = fftw_alloc_complex(N);
	fftw_complex* out = fftw_alloc_complex(N);

	for (int i = 0; i < N; ++i)
	{
		in[i][0] = values[i].real();
		in[i][1] = values[i].imag();;
	}

	fftw_plan plan = fftw_plan_dft_1d(N, in, out, FFTW_BACKWARD,  FFTW_ESTIMATE);
	fftw_execute(plan);

	for (int i = 0; i < N; ++i)
		results[i] = out[i][0] / N;

	fftw_destroy_plan(plan);
	fftw_free(in);
	fftw_free(out);
	fftw_cleanup();

	return results;
}

