#include <iostream>
#include "VectorMath.h"
#include "SF.h"
#include "matplot/matplot.h"

void plotResults(std::vector<std::vector<double>> Matrix)
{
	double trueValue = 1.311028777146120;
	
	std::vector<int> n(Matrix[0].size());
	std::vector<int> sn(Matrix[1].size());
	for (int i = 1; i <= Matrix[0].size(); ++i)
	{
		n[i] = i;
	}

	for (int i = 1; i <= Matrix[1].size(); ++i)
	{
		sn[i] = i;
	}

	Matrix[0] = abs(Matrix[0] - trueValue);
	Matrix[1] = abs(Matrix[1] - trueValue);
	Matrix[2] = abs(Matrix[2] - trueValue);
	
	matplot::semilogy(n, Matrix[0]);
	matplot::hold(matplot::on);
	matplot::semilogy(sn, Matrix[1]);
	matplot::semilogy(n, Matrix[2]);
	matplot::legend({"Trapezoid", "Simpson", "Gaussian Quadrature"});
	matplot::show();

}
