// sample.cpp

#include "PCCCalculator.h"
#include <iostream>

static void PrintResult(const Eigen::MatrixXd &main, const Eigen::MatrixXd &det)
{
	std::cout << "Maintaining directions" << std::endl;
	std::cout << main << std::endl;
	std::cout << "Detatching directions" << std::endl;
	std::cout << det << std::endl;
}

int main(int argc, char argv)
{
	PCCCalculator pcc;
	// simple sample
	Eigen::MatrixXd mat(1, 3);
	mat(0, 0) = 1.0; mat(0, 1) = 0.0; mat(0, 2) = 0.0;
	Eigen::MatrixXd dual_m, dual_d;
	pcc.Dual(dual_m, dual_d, mat);
	PrintResult(dual_m, dual_d);
	Eigen::MatrixXd mat2(0, 3);
	pcc.Dual(dual_m, dual_d, mat, mat2);
	PrintResult(dual_m, dual_d);
	mat2.resize(1, 3);
	mat2(0, 0) = 0.0; mat2(0, 1) = 1.0; mat2(0, 2) = 0.0;
	pcc.Dual(dual_m, dual_d, mat, mat2);
	PrintResult(dual_m, dual_d);
	mat = Eigen::MatrixXd::Identity(3, 3);
	pcc.Dual(dual_m, dual_d, mat);
	PrintResult(dual_m, dual_d);
	return 0;
}