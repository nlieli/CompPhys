#pragma once
#include <vector>
#include <stdexcept>

template <typename T, typename B, typename C>
std::vector<double> linspace(T a, B b, C n = 0)
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
		vector[i] = abs(vector[i]);

	return vector;
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
		vector[i] = scalar - vector[i];

	return vector;
}

template <typename T, typename B>
std::vector<T> vectorScalarSub(B scalar, std::vector<T> vector)
{
	size_t length = vector.size();
	for (size_t i = 0; i < length; ++i)
		vector[i] = scalar - vector[i];

	return vector;
}

template <typename T, typename B>
std::vector<T> operator-(std::vector<T> vector, B scalar)
{
	return vectorScalarSub(vector, scalar);
}

template <typename T, typename B>
std::vector<T> operator-(B scalar, std::vector<T> vector)
{
	return vectorScalarSub(scalar, vector);
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

template <typename T>
std::vector<T> operator<<(std::vector<T> vector, int)
{
	return vectorShiftLeft(vector);
}
