#pragma once
#include <fstream>

enum decimalStyle
{
	comma = ',',
	dot = '.'
};

struct Timer
{
	// properties
	std::chrono::time_point<std::chrono::steady_clock> global_start, global_end;
	std::chrono::duration<float> global_duration;

	std::chrono::time_point<std::chrono::steady_clock> local_start, local_end;
	std::chrono::duration<float> local_duration;
	std::vector<float> laps;

	// methods
	void start();
	void end();
	Timer();
	~Timer();
};

namespace fftw3
{
	std::vector<std::complex<double>> fft(std::vector<double> values, Timer* timer = nullptr);
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


