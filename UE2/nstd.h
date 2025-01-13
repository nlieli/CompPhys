
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
	-scalarProduct - done
	-vectorScalarAdd
	-vectorScalarSub
	-cumsum
	-cumsumVector
	-vectorABS
	-vectorShiftRight
	-vectorShiftLeft
	-arrayErase
-Extend VectorMath to matrices
-create new Matrix specific functions
	-Matrix-Vector Product && Matrix-Matrix Product - done
	-diag iterator type? to find and set values on diags and their offsets - done
-find solution for different type precision
-add rvalue recognition to print function
-print function matrix overload
*/

// std::cerr << "\\033[38;2;255;255;0m[WARNING]: \033[0m" << std::endl;
// std::cerr << "\033[1;31m[ERROR]: \033[0m";

// !!! IMPORTANT !!! Always #include this file last in your project to avoid overload conflicts

#pragma once
#include <type_traits>
#include <iterator>
#include <vector>
#include <optional>
#include <iostream>

/*
You may want to disable the overloads since they are part of the global name space
and thus may interfere with some other functions. To disable set the macro
DISABLE_OVERLOADS to 1
*/

namespace nstd
{
	// custom type traits
	template <typename T>
	struct is_vector : std::false_type {};

	template <typename T>
	struct is_vector<std::vector<T>> : std::true_type {};

	template <typename T>
	struct is_array : std::false_type {};

	template <typename T, size_t N>
	struct is_array<std::array<T, N>> : std::true_type {};

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
	struct is_n_dim_array < std::array<std::array<T, J>, K>> : std::true_type {};


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

	/*
	iterator for n x n matrices that goes over diagonals
	not a standard implementation for an iterator, since
	the container doesnt provide suitable begin() and end()
	functions and i dont want to create a wrapper, they are
	integrated into this iterator type - the iterator is
	fully usable with standard library algorithms such
	as std::distance and std::for_each, etc.
	*/
	template <typename container>
	class diagonal_iterator
	{
		// properties
	public:
		using inner_container = typename container::value_type;
		using value_type = typename inner_container::value_type;
		using pointer = value_type*;
		using reference = value_type&;
		using iterator_category = std::bidirectional_iterator_tag;
		using difference_type = std::ptrdiff_t;

	private:
		size_t m_matrix_size;
		size_t m_row;
		size_t m_row_init;
		size_t m_row_end;
		size_t m_col;
		size_t m_col_init;
		size_t m_col_end;
		container* m_matrix;
		int m_diagonal_identifier;

		// methods
	public:
		diagonal_iterator(container& matrix, int diagonal_identifier = 0)
			: m_matrix(&matrix), m_diagonal_identifier(diagonal_identifier), m_row(0), m_col(diagonal_identifier)
		{
			m_matrix_size = matrix.size();

			if (m_matrix_size != matrix[0].size())
			{
				std::cerr << "\033[38;2;255;255;0m[WARNING]: Using non square matrices may lead to undefined behaviour\033[0m" << std::endl;
			}

			if (std::abs(diagonal_identifier) > m_matrix_size - 1)
			{
				std::cerr << "\033[1;31m[ERROR]: Matrix diagonal is out of bounds!\033[0m" << std::endl;
				throw std::runtime_error("Out of bounds!");
			}

			if (diagonal_identifier >= 0)
			{
				m_row = 0;
				m_col = diagonal_identifier;
				m_row_end = m_matrix_size - diagonal_identifier;
				m_col_end = m_matrix_size;
			}
			else
			{
				m_row = -diagonal_identifier;
				m_col = 0;
				m_row_end = m_matrix_size;
				m_col_end = m_matrix_size + diagonal_identifier;
			}
			m_col_init = m_col;
			m_row_init = m_row;
		}

		reference operator*() { return (*m_matrix)[m_row][m_col]; }

		pointer operator->() { return &(*m_matrix)[m_row][m_col]; }

		diagonal_iterator& operator++()
		{
			++m_row;
			++m_col;
			return *this;
		}

		diagonal_iterator& operator--()
		{
			--m_row;
			--m_col;
			return *this;
		}

		bool operator==(const diagonal_iterator& other) const
		{
			return m_row == other.m_row && m_col == other.m_col;
		}

		bool operator!=(const diagonal_iterator& other) const
		{
			return m_row != other.m_row || m_col != other.m_col;
		}

		difference_type operator-(const diagonal_iterator& other) const
		{
			return static_cast<difference_type>(m_col - other.m_col);
		}

		diagonal_iterator begin() const
		{
			diagonal_iterator begin_iterator = *this;
			begin_iterator.m_col = m_col_init;
			begin_iterator.m_row = m_row_init;
			return begin_iterator;
		}

		diagonal_iterator end() const
		{
			diagonal_iterator end_iterator = *this;
			end_iterator.m_row = m_row_end;
			end_iterator.m_col = m_col_end;
			return end_iterator;
		}
	};

	// make sure to specify T when calling the function
	// for arrays use: auto arr = nstd::linspace<T, std::array<T, points>>(end, start)
	template <typename T, typename container = std::vector<T>>
    container linspace(T start, T end, std::size_t points = 0) 
    {
        static_assert(is_vector_or_array<container>::value, "Container must be an array or vector!");
        double step;
        container result;

        if constexpr (is_array<container>::value)
        {
			if ((points != 0) && points != std::tuple_size_v<container>)
			{
				std::cerr << "\033[1;31m[ERROR]: Container size must equal number of points!" << std::endl;
				throw std::runtime_error("Contianer and point size don't match.");
			}
			else
			{
				points = std::tuple_size_v<container>;
			}
        }


        if (points == 0)
        {
            step = (end > start) ? 1 : -1;
			points = static_cast<std::size_t>(std::abs(end - start) + 1);
        }
        else
        {
            step = static_cast<double>((end - start) / (points - 1));
        }

        if constexpr (is_vector<container>::value)
        {
            result.resize(points);
        }

        for (size_t i = 0; i < points; ++i)
        {
            result[i] = start + i * step;
        }

        return result;
    }

	// math functions
	template <typename T>
	T absNdimArray(const T& container)
	{
		if constexpr (is_vector_or_array<T>::value)
		{
			T result = container;
			for (size_t i = 0; i < container.size(); ++i)
			{
				if constexpr (is_n_dim_array<T>::value)
				{
					result[i] = absNdimArray(container[i]);
				}
				else
				{
					result[i] = std::abs(container[i]);
				}
			}
			return result;
		}
		else
		{
			return std::abs(container);
		}

	}

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

	template <typename T, typename S>
	T scalarAddNdimArray(const T& container, const S& scalar)
	{
		if constexpr (is_vector_or_array<T>::value)
		{
			T result = container;
			for (size_t i = 0; i < container.size(); ++i)
			{
				if constexpr (is_n_dim_array<T>::value)
				{
					result[i] = scalarAddNdimArray(container[i], scalar);
				}
				else
				{
					result[i] = container[i] + scalar;
				}
			}
			return result;
		}
		else
		{
			return container + scalar;
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
	T scalarSubNdimArray(const T& container, const S& scalar)
	{
		if constexpr (is_vector_or_array<T>::value)
		{
			T result = container;
			for (size_t i = 0; i < container.size(); ++i)
			{
				if constexpr (is_n_dim_array<T>::value)
				{
					result[i] = scalarSubNdimArray(container[i], scalar);
				}
				else
				{
					result[i] = container[i] - scalar;
				}
			}
			return result;
		}
		else
		{
			return container - scalar;
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

	template <typename T, typename S>
	T elemPowNdimArray(const T& container, const S& scalar)
	{
		if constexpr (is_vector_or_array<T>::value)
		{
			T result = container;
			for (size_t i = 0; i < container.size(); ++i)
			{
				if constexpr (is_n_dim_array<T>::value)
				{
					result[i] = elemPowNdimArray(container[i], scalar);
				}
				else
				{
					result[i] = std::pow(container[i], scalar);
				}
			}
			return result;
		}
		else
		{
			return std::pow(container, scalar);
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

	template <typename T>
	double norm1dimArray(const T& container)
	{
		return std::sqrt(scalarProduct1dimArray(container, container));
	}

	template <typename T>
	std::vector<std::vector<T>> matrixMultiplication(const std::vector<std::vector<T>>& con1, const std::vector<std::vector<T>>& con2)
	{
		size_t n1 = con1.size();
		size_t m1 = con1[0].size();
		size_t n2 = con2.size();
		size_t m2 = con2[0].size();

		if (m1 != n2)
		{
			std::cerr << "\033[1;31m[ERROR]: Matrices have incompatible sizes!" << std::endl;
			throw std::runtime_error("Incompatible Sizes");
		}

		std::vector<std::vector<T>> result(m1, std::vector<T>(n2));

		for (size_t i = 0; i < m1; ++i)
		{
			for (size_t k = 0; k < n1; ++k)
			{
				for (size_t j = 0; j < n2; ++j)
				{
					result[i][j] += con1[i][k] * con2[k][j];
				}
			}
		}

		return result;
	}

	template <typename T>
	std::vector<T> matrixVectorMultiplication(const std::vector<std::vector<T>>& con1, const std::vector<T>& con2)
	{
		size_t n1 = con1.size();
		size_t m1 = con1[0].size();
		size_t n2 = con2.size();

		if (m1 != n2)
		{
			std::cerr << "\033[1;31m[ERROR]: Matrices have incompatible sizes!" << std::endl;
			throw std::runtime_error("Incompatible Sizes");
		}

		std::vector<T> result(m1);

		for (size_t i = 0; i < m1; ++i)
		{
			for (size_t j = 0; j < n1; ++j)
			{
				result[i] += con1[i][j] * con2[j];
			}
		}

		return result;
	}

	template <typename T>
	std::vector<std::vector<T>> vectorProdToMatrix(std::vector<T>& con1, std::vector<T>& con2)
	{
		size_t rows = con1.size();
		size_t cols = con2.size();
		std::vector<std::vector<T>> result(rows, std::vector<T>(cols));

		for (size_t i = 0; i < rows; ++i)
		{
			for (size_t j = 0; j < cols; ++j)
			{
				result[i][j] = con1[i] * con2[j];
			}
		}

		return result;
	}

	template <typename T>
	std::vector<T> crossProduct(std::vector<T>& con1, std::vector<T>& con2)
	{
		if (con1.size() != 3 && con2.size() != 3) { throw std::runtime_error("Incorrect Size");  }
		std::vector<T> result(3);
		for (size_t i = 0; i < con1.size(); ++i)
			result[i] = con1[(i + 1) % 3] * con2[(i + 2) % 3] - con1[(i + 2) % 3] * con2[(i + 1) % 3];

		return result;
	}


} // end of namespace

// operator overloads
#ifndef DISABLE_NSTD_OVERLOADS 

template <typename T>
T operator+(const T& con1, const T& con2) { return nstd::addNdimArray(con1, con2); }

template <typename T, typename S>
std::enable_if_t<nstd::is_vector_or_array<T>::value && !std::is_same<T, S>::value, T> operator+(const T& container, const S& scalar)
{
	return nstd::scalarAddNdimArray(container, scalar);
}

template <typename T, typename S>
std::enable_if_t<nstd::is_vector_or_array<T>::value && !std::is_same<T, S>::value, T> operator+(const S& scalar, const T& container)
{
	return nstd::scalarAddNdimArray(container, scalar);
}

template <typename T>
T operator-(const T& con1, const T& con2) { return nstd::subNdimArray(con1, con2); }

template <typename T, typename S>
std::enable_if_t<nstd::is_vector_or_array<T>::value && !std::is_same<T, S>::value, T> operator-(const T& container, const S& scalar)
{
	return nstd::scalarSubNdimArray(container, scalar);
}

template <typename T, typename S>
std::enable_if_t<nstd::is_vector_or_array<T>::value && !std::is_same<T, S>::value, T> operator-(const S& scalar, const T& container)
{
	return nstd::scalarSubNdimArray(container, scalar);
}

template <typename T, typename S>
std::enable_if_t<nstd::is_vector_or_array<T>::value && !std::is_same<T, S>::value, T> operator*(const T& container, const S& scalar)
{
	return nstd::scalarMultNdimArray(container, scalar);
}

template <typename T, typename S>
std::enable_if_t<nstd::is_vector_or_array<T>::value && !std::is_same<T, S>::value, T> operator*(const S& scalar, const T& container)
{
	return nstd::scalarMultNdimArray(container, scalar);
}

template <typename T>
std::vector<std::vector<T>> operator*(const std::vector<std::vector<T>>& con1, const std::vector<std::vector<T>>& con2)
{
	return nstd::matrixMultiplication(con1, con2);
}

template <typename T>
std::vector<T> operator*(const std::vector<std::vector<T>>& con1, const std::vector<T>& con2)
{
	return nstd::matrixVectorMultiplication(con1, con2);
}

template <typename T, typename S>
T operator/(const T& container, const S& scalar) { return nstd::scalarDivNdimArray(container, scalar); }

template <typename T>
std::enable_if_t<nstd::is_vector_or_array<T>::value && !nstd::is_n_dim_array<T>::value, double> operator*(const T& con1, const T& con2)
{
	return nstd::scalarProduct1dimArray(con1, con2);
}

template <typename T, typename S>
T operator^(const T& container, const S& scalar)
{
	return nstd::elemPowNdimArray(container, scalar);
}

#endif



