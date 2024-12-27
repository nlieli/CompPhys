/*
TO DO:
-is_vector_or_array  - check
-vector2array and vice versa - later?
-generic print function - check
-VectorMath functions
-Extend VectorMath to matrices
-add diag and other functions
-find solution for different type precision
-add rvalue recognition to print function
-print function matrix overload
*/


#pragma once
#include <type_traits>
#include <vector>

namespace nstd
{
	// custom type traits
	template <typename T>
	struct is_vector_or_array { static constexpr bool value = false; };

	template <typename T>
	struct is_vector_or_array<std::vector<T>> { static constexpr bool value = true; };

	template <typename T, std::size_t N>
	struct is_vector_or_array<std::array<T, N>> { static constexpr bool value = true; };

	template <typename T>
	struct is_n_dim_array { static constexpr bool value = false; };

	template <typename T>
	struct is_n_dim_array<std::vector<std::vector<T>>> { static constexpr bool value = true; };

	template <typename T, size_t J, size_t K>
	struct is_n_dim_array < std::array<std::array<T, J>, K>> { static constexpr bool value = true; };

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
		T result = con1;
		if constexpr (is_n_dim_array<T>::value)
		{
			size_t N = con1.size();
			for (size_t i = 0; i < N; ++i)
				result[i] = addNdimArray(con1[i], con2[i]);
		}
		else if constexpr (is_vector_or_array<T>::value)
		{
			size_t N = con1.size();
			for (std::size_t i = 0; i < N; ++i)
				result[i] = con1[i] + con2[i];
		}
		else
		{
			result = con1 + con2;
		}

		return result;
	}
}

// operator overloads
template <typename T>
T operator+(T con1, T& con2)
{
	return nstd::addNdimArray(con1, con2);
}

