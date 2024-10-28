#pragma once
#include <vector>
#include <cmath>
#include <stdexcept>

// all functions could and should be generalized to iterable types (no time for this atm)

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

template <typename T>
std::vector<T> abs(std::vector<T> vector)
{
	size_t length = vector.size();
	for (size_t i = 0; i < length; ++i)
		vector[i] = std::abs(vector[i]);

	return vector;
}

template <typename T, int size>
std::array<T, size> abs(std::array<T, size> array)
{
	for (T& value : array)
		value = std::abs(value);

	return array;
}


template <typename T, typename B>
std::vector<T> vectorScalarAdd(std::vector<T> vector, B scalar)
{
	size_t length = vector.size();
	for (size_t i = 0; i < length; ++i)
		vector[i] += scalar;

	return vector;
}

template <typename T, typename B>
std::vector<T> vectorScalarAdd(B scalar, std::vector<T> vector)
{
	size_t length = vector.size();
	for (size_t i = 0; i < length; ++i)
		vector[i] += scalar;

	return vector;
}

template <typename T, typename B>
std::vector<T> operator+(std::vector<T> vector, B scalar)
{
	return vectorScalarAdd(vector, scalar);
}

template <typename T, typename B>
std::vector<T> operator+(B scalar, std::vector<T> vector)
{
	return vectorScalarAdd(scalar, vector);
}

template <typename T, typename B>
std::vector<T> vectorScalarSub(std::vector<T> vector, B scalar)
{
	size_t length = vector.size();
	for (size_t i = 0; i < length; ++i)
		vector[i] -= scalar;

	return vector;
}

template <typename T, int size, typename B>
std::array<T, size> arrayScalarSub(std::array<T, size> array, B scalar)
{
	for (T& value : array)
		value -= scalar;

	return array;
}

template <typename T, typename B>
std::vector<T> vectorScalarSub(B scalar, std::vector<T> vector)
{
	size_t length = vector.size();
	for (size_t i = 0; i < length; ++i)
		vector[i] = scalar - vector[i];

	return vector;
}

template <typename T, int size, typename B>
std::array<T, size> scalarArraySub(B scalar, std::array<T, size> array)
{
	for (T& value : array)
		value = scalar - value;

	return array;
}

template <typename T, typename B>
std::vector<T> operator-(std::vector<T> vector, B scalar)
{
	return vectorScalarSub(vector, scalar);
}

template <typename T, int size, typename B>
std::array<T, size> operator-(std::array<T, size> array, B scalar)
{
	return arrayScalarSub(array, scalar);
}

template <typename T, typename B>
std::vector<T> operator-(B scalar, std::vector<T> vector)
{
	return vectorScalarSub(scalar, vector);
}

template <typename T, int size, typename B>
std::array<T, size> operator-(B scalar, std::array<T, size> array)
{
	return scalarArraySub(scalar, array);
}

template <typename T>
std::vector<T> vectorAdd(std::vector<T> v1, std::vector<T> v2) // Addition of two vectors
{
	if (v1.size() != v2.size()) throw std::invalid_argument("Vectors must be the same length");
	size_t length = v1.size();
	for (size_t i = 0; i < length; ++i)
		v1[i] += v2[i];
	return v1;
}

template <typename T>
std::vector<T> operator+(std::vector<T> v1, std::vector<T> v2) // Overload for vector addition	
{
	return vectorAdd(v1, v2);
}

template <typename T>
std::vector<T> vectorSub(std::vector<T> v1, std::vector<T> v2) // Subtraction of two vectors follows common order 
{
	if (v1.size() != v2.size()) throw std::invalid_argument("Vectors must be the same length");
	size_t length = v1.size();
	for (size_t i = 0; i < length; ++i)
		v1[i] -= v2[i];
	return v1;
}

template <typename T>
std::vector<T> operator-(std::vector<T> v1, std::vector<T> v2) // Overload for vector subtraction
{
	return vectorSub(v1, v2);
}

template <typename T, typename B>
std::vector<T> vectorMultiply(std::vector<T> vector, B scalar)
{
	size_t length = vector.size();
	for (size_t i = 0; i < length; ++i)
		vector[i] *= scalar;

	return vector;
}

template <typename T, typename B>
std::vector<T> vectorMultiply(B scalar, std::vector<T> vector)
{
	size_t length = vector.size();
	for (size_t i = 0; i < length; ++i)
		vector[i] *= scalar;

	return vector;
}

template <typename T, typename B>
std::vector<T> operator*(std::vector<T> vector, B scalar)
{
	return vectorMultiply(vector, scalar);
}

template <typename T, typename B>
std::vector<T> operator*(B scalar, std::vector<T> vector)
{
	return vectorMultiply(scalar, vector);
}

template <typename T, typename B>
std::vector<T> vectorDiv(std::vector<T> vector, B scalar)
{
	size_t length = vector.size();
	for (size_t i = 0; i < length; ++i)
		vector[i] /= scalar;
	
	return vector;
}

template <typename T, typename B>
std::vector<T> operator/(std::vector<T> vector, B scalar)
{
	return vectorDiv(vector, scalar);
}

template <typename T, int size>
std::array<T, size> arrayLog10(std::array<T, size> array)
{
	for (T& value : array)
		value = std::log10(value);
	return array;
}

template <typename T>
std::vector<T> vectorShiftRight(std::vector<T> vector)
{
	size_t length = vector.size();
	if (length == 0) return vector;

	for (size_t i = length - 1; i > 0; --i)
		vector[i] = vector[i - 1];

	vector[0] = 0;
	return vector;
}

template <typename T, int size>
std::array<T, size> arrayShiftRight(std::array<T, size> array)
{
	if (size == 0) return array;

	for (int i = size - 1; i > 0; --i)
		array[i] = array[i - 1];

	array[0] = 0;
	return array;
}

template <typename T>
std::vector<T> operator>>(std::vector<T> vector, int)
{
	return vectorShiftRight(vector);
}

template <typename T>
std::vector<T> vectorShiftLeft(std::vector<T> vector)
{
	size_t length = vector.size();
	for (size_t i = 1; i < length; ++i)
		vector[i - 1] = vector[i];

	vector[length - 1] = 0;
}

template <typename T, int size>
std::array<T, size> arrayShiftLeft(std::array<T, size> array)
{
	for (int i = 1; i < size; ++i)
		array[i - 1] = array[i];

	array[size - 1] = 0;
	return array;
}

template <typename T>
std::vector<T> operator<<(std::vector<T> vector, int)
{
	return vectorShiftLeft(vector);
}
