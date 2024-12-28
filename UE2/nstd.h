/*
TO DO:
-is_vector_or_array  - done
-vector2array and vice versa - later?
-generic print function - done
-VectorMath functions
	-vectorAdd - done
	-vectorSub - done
	-vectorScalarMultiply - done
	-vectorScalarDiv - done
	-scalarProduct
	-vectorScalarAdd
	-vectorScalarSub
	-cumsum
	-cumsumVector
	-vectorABS
	-vectorShiftRight
	-vectorShiftLeft
	-arrayErase
-Extend VectorMath to matrices
-add diag and other functions
-find solution for different type precision
-add rvalue recognition to print function
-print function matrix overload
*/


#pragma once
#include <type_traits>
#include <vector>

/*
You may want to disable the overloads since they are part of the global name space
and thus may interfere with some other functions. To disable set the macro 
DISABLE_OVERLOADS to 1
*/
#define DISABLE_OVERLOADS 0 

namespace nstd
{
	// custom type traits
	template <typename T>
	struct is_vector_or_array : std::false_type {};

	template <typename T>
	struct is_vector_or_array<std::vector<T>> : std::true_type {};

	template <typename T, std::size_t N>
	struct is_vector_or_array<std::array<T, N>> : std::true_type {};

	template <typename T>
	struct is_n_dim_array : std::false_type {};

	template <typename T>
	struct is_n_dim_array<std::vector<std::vector<T>>> : std::true_type {};

	template <typename T, size_t J, size_t K>
	struct is_n_dim_array < std::array<std::array<T, J>, K>> :std::true_type {};

	// utility functions
	template <typename T>
	void print(const T& container, const char* variableName = nullptr)
	{
		if constexpr (is_vector_or_array<T>::value)
		{
			size_t size = container.size();

			if (variableName)
				std::cout << variableName << " = ";
			std::cout << "[";

			if constexpr (is_n_dim_array<T>::value) // maybe stupid and unnecessary for a print function (proof of concept for other operations)
			{
				for (size_t i = 0; i < size; ++i)
					nstd::print(container[i]);
			}
			else
			{
				for (size_t i = 0; i < size - 1; ++i)
					std::cout << container[i] << ", ";
				std::cout << container[size - 1]
					<< "]" << std::endl;
			}
		}
		else
		{
			if (variableName)
				std::cout << variableName << " = ";
			std::cout << container << std::endl;
		}
	}

	// math functions
	template <typename T>
	T addNdimArray(const T& con1, const T& con2)
	{
		if constexpr (is_vector_or_array<T>::value)
		{
			if (con1.size() != con2.size()) {
				std::cerr << "\033[1;31m[ERROR] Containers must be the same size!\033[0m";
				throw std::runtime_error("Incompatible sizes!");
			}

			T result = con1;
			for (size_t i = 0; i < con1.size(); ++i)
			{
				if constexpr (is_n_dim_array<T>::value)
				{
					result[i] = addNdimArray(con1[i], con2[i]);
				}
				else
				{
					result[i] = con1[i] + con2[i];
				}
			}
			return result;
		}
		else
		{
			return con1 + con2;
		}
	}

	template <typename T>
	T subNdimArray(const T& con1, const T& con2)
	{
		if constexpr (is_vector_or_array<T>::value)
		{
			if (con1.size() != con2.size()) {
				std::cerr << "\033[1;31m[ERROR] Containers must be the same size!\033[0m";
				throw std::runtime_error("Incompatible sizes!");
			}

			T result = con1;
			for (size_t i = 0; i < con1.size(); ++i)
			{
				if constexpr (is_n_dim_array<T>::value)
				{
					result[i] = subNdimArray(con1[i], con2[i]);
				}
				else
				{
					result[i] = con1[i] - con2[i];
				}
			}
			return result;
		}
		else
		{
			return con1 - con2;
		}
	}

	template <typename T, typename S>
	T scalarMultNdimArray(const T& container, const S& scalar) 
	{
		if constexpr (is_vector_or_array<T>::value)
		{
			T result = container;
			for (size_t i = 0; i < container.size(); ++i)
			{
				if constexpr (is_n_dim_array<T>::value)
				{
					result[i] = scalarMultNdimArray(container[i], scalar);
				}
				else
				{
					result[i] = container[i] * scalar;
				}
			}
			return result;
		}
		else
		{
			return container * scalar;
		}
	}

	template <typename T, typename S>
	T scalarDivNdimArray(const T& container, const S& scalar)
	{
		if constexpr (is_vector_or_array<T>::value)
		{
			T result = container;
			for (size_t i = 0; i < container.size(); ++i)
			{
				if constexpr (is_n_dim_array<T>::value)
				{
					result[i] = scalarDivNdimArray(container[i], scalar);
				}
				else
				{
					result[i] = container[i] / scalar;
				}
			}
			return result;
		}
		else
		{
			return container / scalar;
		}
	}

	template <typename T>
	std::enable_if_t<is_vector_or_array<T>::value && !is_n_dim_array<T>::value, double> scalarProduct1dimArray(const T& con1, const T& con2) // also works on matrix hmmmm...
	{
			size_t N = con1.size();
			if (N != con2.size()) {
				std::cerr << "\033[1;31m[ERROR] Containers must be the same size!\033[0m";
				throw std::runtime_error("Incompatible sizes!");
			}

			double result = 0;
			for (size_t i = 0; i < N; ++i)
				result += con1[i] * con2[i];

			return result;
	}

} // end of namespace

// operator overloads
#if DISABLE_OVERLOADS == 0
template <typename T>
T operator+(const T& con1, const T& con2) { return nstd::addNdimArray(con1, con2); }

template <typename T>
T operator-(const T& con1, const T& con2) { return nstd::subNdimArray(con1, con2); }

template <typename T, typename S>
std::enable_if_t<nstd::is_vector_or_array<T>::value && !std::is_same<T, S>::value, T> operator*(const T& container, const S& scalar) 
{ return nstd::scalarMultNdimArray(container, scalar); }

template <typename T, typename S>
std::enable_if_t<nstd::is_vector_or_array<T>::value && !std::is_same<T, S>::value, T> operator*(const S& scalar, const T& container) 
{ return nstd::scalarMultNdimArray(container, scalar); }

template <typename T, typename S>
T operator/(const T& container, const S& scalar) { return nstd::scalarDivNdimArray(container, scalar); }

template <typename T>
std::enable_if_t<nstd::is_vector_or_array<T>::value && !nstd::is_n_dim_array<T>::value, double> operator*(const T& con1, const T& con2) 
{ return nstd::scalarProduct1dimArray(con1, con2); }

#endif



