#include <iostream>
#include "SF.h"
#include "matplot/matplot.h"

void plotResults(std::vector<std::vector<double>> Matrix)
{
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

	printVector(sn);
	printVector(Matrix[2]);

	//matplot::plot(n, Matrix[0]);
	//matplot::hold(matplot::on);
	//matplot::plot(sn, Matrix[1]);
	//matplot::plot(n, Matrix[2]);
	//matplot::show();

}
