#pragma once
#include <vector>
#include <iostream>

unsigned int factorial(unsigned int x);

std::vector<double>& nthlegendre(size_t n);

template <typename T>
void printVector(std::vector<T>& vector);
