#pragma once
#include <vector>
#include <array>
#include "SF.h"

namespace SF{
	class bracketInstance
	{
	private:
		std::vector<double> m_function;
		double m_split;
		double epsilon = 1e-4;
		unsigned int m_functionOrder;

	public:
		double m_upperBoundary;
		double m_lowerBoundary;
		bool m_rootState;
		double m_intervalLength;
		unsigned int m_Count;
		bracketInstance* m_parent;
		bracketInstance* m_SIP1;
		bracketInstance* m_SIP2;


	private:
		bool checkBracket();

		void createSubinstance()
		{
			m_split = (m_lowerBoundary + m_upperBoundary) / 2;
			bracketInstance subinstance1(m_function, m_lowerBoundary, m_split, this);
			bracketInstance subinstance2(m_function, m_split, m_lowerBoundary, this);
			m_SIP1 = &subinstance1;
			m_SIP2 = &subinstance2;
		}

		void updateCount()
		{

		}


	public:
		bracketInstance(std::vector<double> function, double a, double b, bracketInstance* parent = nullptr)
		{
			m_upperBoundary = b;
			m_lowerBoundary = a;
			m_function = function;
			m_parent = parent;
			m_intervalLength = abs(b - a);
			m_Count = 0;
			m_functionOrder = function.size() - 1;
			m_rootState = checkBracket();
			m_SIP1 = nullptr;
			m_SIP2 = nullptr;
		}

		void subinstanceUpdate()
		{

		}
	};
}
